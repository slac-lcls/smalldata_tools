import io
import os
import re

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
      <a href="../{prevPage}/report.html" style="margin-right:150px"> Previous </a>

      <a href="../{nextPage}/report.html"> Next </a>
    </div>
"""


def getFormattedPageTitle(experiment: str, run: int, pageTitle: str) -> str:
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
    prevPage: str = getFormattedPageTitle(experiment, run - 1, pageTitle)
    nextPage: str = getFormattedPageTitle(experiment, run + 1, pageTitle)
    if "/" in prevPage:
        # HANDLE pageTitle formats of type: Summary/Summary_RUN
        prevPage = prevPage.split("/")[-1]
        nextPage = nextPage.split("/")[-1]

    return LINK_DIV.format(prevPage=prevPage, nextPage=nextPage)


def prepareHtmlReport(tabs: pn.Tabs, experiment: str, run: int, pageTitle: str) -> None:
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
