import math
import os
import time
import datetime
import urllib
try:
    from urllib.request import urlopen
except:
    from urllib2 import urlopen
import matplotlib.pyplot as plt
import json

base_url = "http://pscaa02"
retrieval_url = base_url + ":17668/retrieval/data/getData.json"
mgmt_url = base_url + ":17665/mgmt/bpl/"
pv_arg = "?pv={}"
url_arg = "&{0}={1}"
url_flag = "&{0}"
date_spec_format   = "{0:04}-{1:02}-{2:02}T{3:02}:{4:02}:{5:02}.{6:03}Z"

class EpicsArchive(object):
    """
    Class that accesses data from the new archiver.
    Currently supports getting points, plotting points, and searching for pvs.
    """
    def __init__(self):
        self._pts_cache = None
        self._pv_cache = None

    def get_points(self, PV=None, start=30, end=None, unit="days", chunk=False, two_lists=False, raw=False):
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
                json_obj = self._get_json(PV, json_start, json_end, chunk)
                pts = self._json_to_pts(json_obj)
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

    def plot_points(self, PV=None, start=30, end=None, unit="days", chunk=False):
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
        plt.step(t, x, where='post', color='k')
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
            There are three ways to provide these time periods:
            1. numeric arguments are relative time: e.g. start=30, unit="days"
                begins the archiver request at 30 days ago.
            2. list/tuple arguments are absolute time, in the format
                [year, month, day, hour, min, sec], which does not need to be
                fully specified. For example, start=[2016,1,14] end=[2016,2,14]
                gives you the one month period between Jan 14 and Feb 14.
            3. None gives you the endpoint. start=None is the beginning of
                archivable data, end=None is now.

        The unit argument is only used if start or end is an integer. Most
            units of time are valid, expressed as a string, and many aliases
            are supported.
        """
        pass # This only exists to hold the common doc string elements
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

        # Convert to datetime
        start = to_datetime(start, unit)
        end = to_datetime(end, unit)

        # Adjust to UTC
        utcnow = datetime.datetime.utcnow()
        herenow = datetime.datetime.now()
        offset = round((utcnow - herenow).total_seconds()/(60*60), 2)

        utc_start = start + datetime.timedelta(hours=offset)
        utc_end = end + datetime.timedelta(hours=offset)

        # Switch to date arrays for url construction
        json_start = datetime_to_array(utc_start)
        json_end = datetime_to_array(utc_end)

        return json_start, json_end

    def _get_json(self, PV, start, end, chunk):
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
        url += pv_arg.format(urllib.parse.quote(PV, safe=""))
        url += url_arg.format("from", date_format(*start))
        url += url_arg.format("to", date_format(*end))
        if not chunk:
            url += url_flag.format("donotchunk")
        return url_query(url)

    def _json_to_pts(self, json_obj):
        """
        Inteprets a data retrieval json object as an array of tuple points.
        """
        return [ (x["secs"], x["val"]) for x in json_obj[0]["data"] ]

    def _pts_to_arrays(self, pts):
        """
        Changes an array of tuples into two lists.
        """
        t_array = []
        val_array = []
        for (t, val) in pts:
            t_array.append(t)
            val_array.append(val)
        return t_array, val_array

days_map = {}
days_map.update({ x: 365              for x in ("years", "year", "yr", "y")        })
days_map.update({ x: 365./12          for x in ("months", "month", "mon")          })
days_map.update({ x: 7                for x in ("weeks", "week", "wks", "wk", "w") })
days_map.update({ x: 1                for x in ("days", "day", "d")                })
days_map.update({ x: 1./24            for x in ("hours", "hour", "hrs", "hr", "h") })
days_map.update({ x: 1./24/60         for x in ("minutes", "mins", "min")          })
days_map.update({ x: 1./24/60/60      for x in ("seconds", "secs", "sec", "s")     })
days_map.update({ x: 1./24/60/60/1000 for x in ("milliseconds", "msec", "ms")      })

def url_query(url):
    """Makes the URL request."""
    req = urlopen(url)
    data = json.load(req)
    return data

def to_datetime(arg, unit):
    """
    Convert some form of date input into a datetime object.
    """
    if isinstance(arg, datetime.datetime):
        return arg
    if isinstance(arg, (int, float)):
        if (arg < 10000):
            return datetime_ago(arg, unit)
        else: #big - must be time in secs since epoch.
            return datetime.datetime.fromtimestamp(int(arg))
    if isinstance(arg, (list, tuple)):
        arg = list(arg)
        while len(arg) < 3:
            arg.append(0)
        return datetime.datetime(*arg)
 
def datetime_ago(delta, unit):
    """
    Return datetime object at the time of delta units ago.

    For example, datetime_ago(2, "days") would return a datetime object at the
    time equivalent to 2 days ago.
    """
    days = delta * days_map[unit]
    return datetime.datetime.now() - datetime.timedelta(days)

def datetime_to_array(dt):
    """
    Convert datetime object to an array that can be passed to date_format.
    """
    return [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]

def date_format(year=2015, month=1, day=1, hr=0, min=0, s=0, ms=0):
    """Convert date/time parameters to date format string for archiver"""
    d = date_spec_format.format(year, month, day, hr, min, s, ms)
    return urllib.parse.quote(d, safe="")

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
    _, term_width = os.popen('stty size', 'r').read().split()
    n_cols = int(term_width)/col_width
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

