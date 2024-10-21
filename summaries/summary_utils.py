"""Includes functions for creating summary plots.

Functions
---------
getFormattedPageTitle(experiment: str, run: int, pageTitle: str) -> str: Format a
    an uncompiled summary page title to include experiment/run etc.

getLinkDiv(experiment: str, run: int, pageTitle: str) -> str: Update the navigation
    link portion (LINK_DIV) with the correct links.

updateSymLinks(fullPath: str, run: int) -> None: Update the symbolic links pointing
    to the most recent and first run summary report HTML.

prepareHtmlReport(tabs: pn.Tabs, experiment: str, run: int, pageTitle: str) -> None:
    Prepare the full HTML report including navigation, update symbolic links and write
    the report file.
"""

import io
import os
import re
from typing import Optional

import panel as pn


LINK_STYLE: str = """
    <style>
      a:link, a:visited {
        background-color: #FF0033;
        border-radius: 8px;
        border-style: none;
        box-sizing: border-box;
        color: #FFFFFF;
        cursor: pointer;
        display: inline-block;
        font-family: "Haas Grot Text R Web", "Helvetica Neue", Helvetica, Arial, sans-serif;
        font-size: 14px;
        font-weight: 500;
        height: 40px;
        line-height: 20px;
        list-style: none;
        margin: 0;
        outline: none;
        padding: 10px 16px;
        position: relative;
        text-align: center;
        text-decoration: none;
        transition: color 100ms;
        vertical-align: baseline;
        user-select: none;
        -webkit-user-select: none;
        touch-action: manipulation;
      }
      a:hover, a:active {
        background-color: #F082AC;
      }
    </style>
"""

LINK_DIV: str = """
    <div style:"width: 100%; border:1px solid red">
      <a href="{currPage}/../.first.html" style="margin-right:150px"> First </a>

      <a href="{prevPage}/report.html" style="margin-right:150px"> Previous </a>

      <a href="{nextPage}/report.html" style="margin-right:150px"> Next </a>

      <a href="{currPage}/../.latest.html"> Latest </a>
    </div>
    <br><br>
"""


def getFormattedPageTitle(experiment: str, run: int, pageTitle: str) -> str:
    """Apply formatting to a run summary page title.

    Parameters
    ----------
    experiment (str) Name of the experiment.

    run (int) Current run.

    pageTitle (str) Summary page title with format specifiers un-compiled.
        E.g. BeamlineSummary/BeamlineSummary_Run{run:04d}.

    Returns
    -------
    formattTitle (str) Page title with formatting applied.
    """
    runPtn: str = "{run.*}"
    formattedTitle: str = re.sub(
        runPtn, lambda match: match[0].format(run=run), pageTitle
    )

    hutch: str = experiment[:3]
    hutchPtn: str = "{hutch}"
    formattedTitle = re.sub(
        hutchPtn, lambda match: match[0].format(hutch=hutch), formattedTitle
    )

    expmtPtn: str = "{experiment}"
    formattedTitle = re.sub(
        expmtPtn, lambda match: match[0].format(experiment=experiment), formattedTitle
    )

    return formattedTitle


def getLinkDiv(experiment: str, run: int, pageTitle: str) -> str:
    """For an experiment, run and summary page title, update button links.

    Parameters
    ----------
    experiment (str) Experiment name.

    run (int) Run number.

    pageTitle (str) Summary page title with format specifiers un-compiled.
        E.g. BeamlineSummary/BeamlineSummary_Run{run:04d}.

    Returns
    -------
    link_div (str) Formatted HTML for the <div> containing navigation links.
    """
    hutch: str = experiment[:3].upper()
    URLBase: str = "https://pswww.slac.stanford.edu/experiment_results/"
    expBase: str = f"{URLBase}/{hutch}/{experiment}"

    prevPageTitle: str = getFormattedPageTitle(experiment, run - 1, pageTitle)
    currPageTitle: str = getFormattedPageTitle(experiment, run, pageTitle)
    nextPageTitle: str = getFormattedPageTitle(experiment, run + 1, pageTitle)
    prevPage: str = f"{expBase}/{prevPageTitle}"
    currPage: str = f"{expBase}/{currPageTitle}"
    nextPage: str = f"{expBase}/{nextPageTitle}"

    return LINK_DIV.format(prevPage=prevPage, currPage=currPage, nextPage=nextPage)


def updateSymLinks(fullPath: str, run: int) -> None:
    """Update symbolic links pointing to most recent and first run summary.

    WARNING: This function only works for summaries that are alone in a parent
    folder. Either they are nested, for example, like Summary/Summary_Run0000,
    OR, they are at the top level of the experiment's stats/summary directory
    but no other types of summary plots have been created.

    This function also assumes the ONLY NUMBER in the page titles is the RUN
    NUMBER.

    Parameters
    ----------
    fullPath (str) Full path to the CURRENT run summary being processed.

    run (int) The current run being processed.
    """
    currentSummaries: list[str] = [
        d for d in os.listdir(f"{fullPath}/..") if d[0] != "."
    ]
    currentRuns: list[int] = []
    for summ in currentSummaries:
        res: Optional[re.Match] = re.search(r"[0-9]+", summ)
        if res is not None:
            currentRuns.append(int(res[0]))

    currentRuns.sort()
    if len(currentRuns) > 0:
        if run <= currentRuns[0]:
            if os.path.exists(f"{fullPath}/../.first.html"):
                os.remove(f"{fullPath}/../.first.html")
            os.symlink(f"{fullPath}/report.html", f"{fullPath}/../.first.html")
        if run >= currentRuns[-1]:
            if os.path.exists(f"{fullPath}/../.latest.html"):
                os.remove(f"{fullPath}/../.latest.html")
            os.symlink(f"{fullPath}/report.html", f"{fullPath}/../.latest.html")


def prepareHtmlReport(tabs: pn.Tabs, experiment: str, run: int, pageTitle: str) -> None:
    """Given run summary plots prepare an HTML report including navigation links.

    Parameters
    ----------
    tabs (pn.Tabs) Tabular run summary plots to include in HTML report.

    experiment (str) Experiment name.

    run (int) Run number.

    pageTitle (str) Summary page title with format specifiers un-compiled.
        E.g. BeamlineSummary/BeamlineSummary_Run{run:04d}.
    """
    f: io.BytesIO = io.BytesIO()
    tabs.save(f)
    f.seek(0)
    tabsHtmlBytes: bytes = f.read()

    htmlOut: str = ""
    for line in tabsHtmlBytes.decode("UTF-8").split("\n"):
        htmlOut += f"{line}\n"
        if "<title>" in line:
            htmlOut += LINK_STYLE
        elif "<body>" in line:
            htmlOut += getLinkDiv(experiment, run, pageTitle)

    hutch: str = experiment[:3]
    elogBaseDir: str = f"/sdf/data/lcls/ds/{hutch}/{experiment}/stats/summary"
    fullPath: str = f"{elogBaseDir}/{getFormattedPageTitle(experiment, run, pageTitle)}"

    if not os.path.isdir(fullPath):
        os.makedirs(fullPath)
        print("Made Directory to save data:", fullPath)

    with open(f"{fullPath}/report.html", "w") as report:
        report.write(htmlOut)

    # Call updateSymLinks after making new summary page so calculation works
    updateSymLinks(fullPath, run)
