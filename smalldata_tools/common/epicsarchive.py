import math
import os
import time
import datetime
from dateutil.relativedelta import relativedelta
from dateutil import tz
import urllib
import asyncio
import httpx
import matplotlib.pyplot as plt
import json
import logging

logger = logging.getLogger(__name__)


base_url = "https://pswww.slac.stanford.edu/archiveviewer"
retrieval_url = base_url + "/retrieval/data/getData.json"
mgmt_url = base_url + ":17665/mgmt/bpl/"
pv_arg = "?pv={}"
# url_arg = "\&{0}={1}"
url_arg = "&{0}={1}"
url_flag = "&{0}"
date_spec_format = "{0:04}-{1:02}-{2:02}T{3:02}:{4:02}:{5:02}.{6:03}Z"


class EpicsArchive(object):
    """
    Class that accesses data from the new archiver.
    Currently supports getting points, plotting points, and searching for pvs.
    """

    def __init__(self):
        self._pts_cache = None
        self._pv_cache = None

    def get_point(self, PV=None, when=None, value_only=False):
        """
        Get the value from a particular point in time. Returns (timestamp, value)
        unless value_only is true.  Returns None if the archiver does not have
        a value for that time.

        when should either be a timestamp or datetime object.

        This tries to be a bit smarter about detecting disconnects.
        """
        when_start = to_datetime(when, "days")
        json_start, json_end = self._json_args(when, when, "days")
        incr = datetime.timedelta(
            seconds=60
        )  # How far to look forward for a second data point.
        found = False
        now = datetime.datetime.now().astimezone(datetime.timezone.utc)
        while not found:
            when_end = when_start + incr
            if when_end > now:
                when_end = now
            json_end = datetime_to_array(when_end)
            json_obj = self._get_json(PV, json_start, json_end, False)
            if len(json_obj) == 0:  # Time must be before archiving!
                return None
            if len(json_obj[0]["data"]) > 1 or when_end >= now:
                found = True
            else:
                # We only have 0-1 data points.
                if incr.days == 0:
                    incr = 10 * incr
                else:
                    incr = incr + datetime.timedelta(days=1)
        data = json_obj[0]["data"]
        # Now, data[0] is a data point archived before when. But
        # data[1] might contain a fields/cnxlostepsecs tag that
        # invalidates it!
        #
        # What if we *don't* have a second point though.  Does this
        # mean that we are *still* disconnected?  Or is the data valid?
        #
        # For now, let's pretend that it's valid.
        if (
            len(data) > 1
            and "fields" in data[1].keys()
            and "cnxlostepsecs" in data[1]["fields"].keys()
            and float(data[1]["fields"]["cnxlostepsecs"]) < when_start.timestamp()
        ):
            return None
        if value_only:
            return data[0]["val"]
        else:
            return (data[0]["secs"] + data[0]["nanos"] / 1e9, data[0]["val"])

    async def get_points(
        self,
        PV=None,
        start=30,
        end=None,
        unit="days",
        chunk=False,
        two_lists=False,
        raw=False,
        useMS=False,
    ):
        """
        Get points from the archive, returning them as a list of tuples.
        You may set two_lists=True to get a list of positions and a list of
            times instead of a single list of tuples.
        You may set raw=True to get times as seconds from the epoch instead of
            as strings.
        """
        if PV is None:
            pts, _ = self._check_cache()
            if pts is None:
                return
        else:
            json_start, json_end = self._json_args(start, end, unit)
            if valid_date_arrays(json_start, json_end):
                json_obj = await self._get_json(PV, json_start, json_end, chunk)
                if isinstance(json_obj, int):  # Handle http error code
                    logger.warning(
                        f"Not able to retrieve PV {PV} (http code {json_obj})."
                    )
                    return []
                pts = self._json_to_pts(json_obj, useMS)
                self._pts_cache = pts
                self._pv_cache = PV
            else:
                print("Invalid dates! Start must be BEFORE end!")
                pts = []
        if not raw:
            pts = pts_string_time(pts)
        if two_lists:
            return self._pts_to_arrays(pts)
        else:
            return pts

    def plot_points(self, PV=None, start=5, end=None, unit="days", chunk=False):
        """
        Use matplotlib to plot points from the archive.
        """
        if PV is None:
            pts, PV = self._check_cache()
            if pts is None:
                return
        else:
            pts = self.get_points(PV, start, end, unit, chunk, False, True)
        t, x = self._pts_to_arrays(pts)
        for i in range(len(t)):
            t[i] = datetime.datetime.fromtimestamp(t[i])
        fig = plt.figure()
        fig.canvas.set_window_title("EPICS Archive " + PV)
        plt.plot(t, x, "k.")
        plt.step(t, x, where="post", color="k")
        plt.ylabel(PV)
        plt.show(block=False)

    def _check_cache(self):
        """
        Return cache if it exists, otherwise print a message.
        """
        if self._pts_cache is None:
            print("No cached points to use. Please choose a PV.")
            return None, None
        else:
            print("Found cached data for {}".format(self._pv_cache))
            return self._pts_cache, self._pv_cache

    def search_pvs(self, glob, do_print=True):
        """
        Queries the archiver with a PV search using glob format.

        If do_print=True, prints the PVs nicely on your screen.
        If do_print=False, returns the list of matches.
        """
        url = mgmt_url + "getAllPVs"
        url += pv_arg.format(urllib.parse.quote(glob, safe=""))
        pvs = url_query(url)
        if do_print:
            success = list_print(pvs)
            if not success:
                print("No PVs found.")
        else:
            return pvs

    def __interface(self):
        """
        Arguments:

        PV is the string PV to look up in the archiver. If PV is None, we'll
            use the points found in the most recent call to get_points or
            plot_points.

        start and end represent the time period to use.
            There are four ways to provide these time periods:
            1. numeric arguments are relative time: e.g. start=30, unit="days"
                begins the archiver request at 30 days ago.
            2. Large integers (>=10000) give time since epoch.
            3. list/tuple arguments are absolute time, in the format
                [year, month, day, hour, min, sec], which does not need to be
                fully specified. For example, start=[2016,1,14] end=[2016,2,14]
                gives you the one month period between Jan 14 and Feb 14.
            4. None gives you the endpoint. start=None is the beginning of
                archivable data, end=None is now.

        The unit argument is only used if start or end is an integer. Most
            units of time are valid, expressed as a string, and many aliases
            are supported.
        """
        pass  # This only exists to hold the common doc string elements

    get_points.__doc__ += __interface.__doc__
    plot_points.__doc__ += __interface.__doc__

    def _json_args(self, start, end, unit):
        """
        Change the user interface arguments into arguments suitable for the
        URL json request.
        """
        # Default endpoints
        if start is None:
            start = [2010]
        if end is None:
            end = datetime_to_array(datetime.datetime.now())

        # Convert to UTC datetime
        start = to_datetime(start, unit)
        end = to_datetime(end, unit)

        # Switch to date arrays for url construction
        json_start = datetime_to_array(start)
        json_end = datetime_to_array(end)

        return json_start, json_end

    async def _get_json(self, PV, start, end, chunk):
        """
        Do a url query of the new archiver. Return the json result unmodified.

        PV: string PV to look up in the archiver.
        start, end: array of the form [year, month, day, hour, min, sec, ms]
            where you may omit later entries as desired. For example, [2015]
            is a valid argument and would be the start of 2015. [2016, 4, 15]
            is a valid argument and would be the start of April 15, 2016.
        chunk: boolean for whether or not you want data to be chunked.
        """
        url = retrieval_url
        url += pv_arg.format(PV)
        url += url_arg.format("from", date_format(*start))
        url += url_arg.format("to", date_format(*end))
        if not chunk:
            url += url_flag.format("donotchunk")
        logger.debug(f"URL: {url}")
        data = await url_query(url)
        return data

    def _json_to_pts(self, json_obj, useMS):
        """
        Interprets a data retrieval json object as an array of tuple points.
        """
        return [
            (x["secs"] + (x["nanos"] / 1e9 if useMS else 0), x["val"])
            for x in json_obj[0]["data"]
        ]

    def _pts_to_arrays(self, pts):
        """
        Changes an array of tuples into two lists.
        """
        t_array = []
        val_array = []
        for t, val in pts:
            t_array.append(t)
            val_array.append(val)
        return t_array, val_array


days_map = {}
days_map.update({x: 365 for x in ("years", "year", "yr", "y")})
days_map.update({x: 365.0 / 12 for x in ("months", "month", "mon")})
days_map.update({x: 7 for x in ("weeks", "week", "wks", "wk", "w")})
days_map.update({x: 1 for x in ("days", "day", "d")})
days_map.update({x: 1.0 / 24 for x in ("hours", "hour", "hrs", "hr", "h")})
days_map.update({x: 1.0 / 24 / 60 for x in ("minutes", "mins", "min")})
days_map.update({x: 1.0 / 24 / 60 / 60 for x in ("seconds", "secs", "sec", "s")})
days_map.update({x: 1.0 / 24 / 60 / 60 / 1000 for x in ("milliseconds", "msec", "ms")})


async def url_query(url):
    """
    Makes the URL request.
    Comment: httpx automatically encodes the URL, so it does not need
    to be done here (with urllib.parse.quote for example).
    https://www.w3schools.com/tags/ref_urlencode.ASP
    """
    async with httpx.AsyncClient() as client:
        req = await client.get(url)
    if req.status_code != 200:
        return req.status_code
    return req.json()


def to_datetime(arg, unit):
    """
    Convert some form of date input into a UTC datetime object.
    """
    if isinstance(arg, datetime.datetime):
        return arg.astimezone(datetime.timezone.utc)
    if isinstance(arg, (int, float)):
        if arg < 10000:
            return datetime_ago(arg, unit)
        else:  # big - must be time in secs since epoch.
            return datetime.datetime.fromtimestamp(int(arg)).astimezone(
                datetime.timezone.utc
            )
    if isinstance(arg, (list, tuple)):
        arg = list(arg)
        while len(arg) < 3:
            arg.append(0)
        return datetime.datetime(*arg).astimezone(datetime.timezone.utc)


def datetime_ago(delta, unit):
    """
    Return datetime object at the time of delta units ago.

    For example, datetime_ago(2, "days") would return a datetime object at the
    time equivalent to 2 days ago.
    """
    days = delta * days_map[unit]
    t = datetime.datetime.now() - datetime.timedelta(days)
    return t.astimezone(datetime.timezone.utc)


def datetime_to_array(dt):
    """
    Convert UTC datetime object to an array that can be passed to date_format.
    """
    return [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second, dt.microsecond]


def date_format(year=2015, month=1, day=1, hr=0, min=0, s=0, ms=0):
    """Convert date/time parameters to date format string for archiver"""
    d = date_spec_format.format(year, month, day, hr, min, s, ms)
    # return urllib.parse.quote(d, safe="")  # if urllib is used
    return d


def valid_date_arrays(start, end):
    """
    Checks if start is temporally before end.
    """
    start_copy, end_copy = list(start), list(end)
    for arg in (start_copy, end_copy):
        while len(arg) < 3:
            arg.append(1)
    time_delta = datetime.datetime(*end_copy) - datetime.datetime(*start_copy)
    if time_delta > datetime.timedelta(0):
        return True
    return False


def pts_string_time(pts):
    """
    Convert array of points of the form (unix timestamp, value) to the form
    (ctime string, value)
    """
    newpts = []
    for i in range(len(pts)):
        newpts.append((time.ctime(pts[i][0]), pts[i][1]))
    return newpts


def list_print(data):
    """
    Prints list data in columns as you'd expect from a terminal.
    """
    if len(data) == 0:
        return False
    text_data = [str(i) for i in data]
    col_width = max([len(i) for i in text_data]) + 2
    _, term_width = os.popen("stty size", "r").read().split()
    n_cols = int(term_width) / col_width
    n_text = len(text_data)
    n_full_col = int(math.ceil(float(n_text) / n_cols))
    text_rows = []
    for i in range(n_full_col):
        text_rows.append([])
    row = 0
    for text in text_data:
        text_rows[row].append(text)
        if row < n_full_col - 1:
            row += 1
        else:
            row = 0
    line_elem = "{0:{1}}"
    for row in text_rows:
        line = ""
        for text in row:
            line += line_elem.format(text, col_width)
        print(line)
    return True


def ts_to_datetime(ts, timezone="UTC"):
    """
    The LCLS timestamp is using the EPICS epoch, which is
    20 years after the POSIX epoch.
    """
    sec = ts >> 32 & 0xFFFFFFFF
    nsec = nsec = ts & 0xFFFFFFFF
    timezone = tz.gettz(timezone)
    dt = datetime.datetime.fromtimestamp(sec + nsec * 1e-9, tz=timezone)
    dt += relativedelta(years=20)
    return dt
