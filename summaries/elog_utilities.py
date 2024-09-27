"""
A series of generic functions for interacting with various eLog API endpoints.

These functions facilitate eLog actions such as posting messages, run tables,
retrieving tag data and others. More specific eLog functions, e.g. to post a
*specific* message should be included elsewhere and constructed from the tools
available here.

Functions:
    getElogBasicAuth(exp) -> HTTPBasicAuth
        Returns a HTTP authentication object for *ACTIVE* experiments assuming
        an "opr" account (xppopr, mfxopr...).
    getRunsWithTag(exp, tag, http_auth) -> List[int]
        Returns a list of runs for a specific experiment with a corresponding
        tag.
    postElogMsg(exp, msg, tag, title, files)
        Post a message to the eLog, optionally with a tag, title, or file
        attachments.
    postRunTable(runtable_data)
        Post data to a Run Table.
"""

__all__ = ["getElogBasicAuth", "getRunsWithTag", "postElogMsg", "postRunTable"]

import os
import logging
import mimetypes
import socket
from typing import Optional, List, Union, Tuple, Dict, Any

import requests
from requests.auth import HTTPBasicAuth

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

BASE_URL: str = "https://pswww.slac.stanford.edu/ws-auth/lgbk/lgbk"

def getKerberosAuthHeaders() -> dict: ...

def getElogBasicAuth(exp: str) -> HTTPBasicAuth:
    """Return an authentication object for the eLog API for an opr account.

    This method will only work for active experiments. "opr" accounts are
    removed from the authorized users list after the experiment ends.

    Paramters
    ---------
    exp (str) Experiment name (to determine operator username).

    Returns
    -------
    http_auth (HTTPBasicAuth) Authentication for eLog API.
    """
    opr_name: str = f"{exp[:3]}opr"
    hostname: str = socket.gethostname()
    if hostname.find('sdf') >= 0:
        auth_path: str = "/sdf/group/lcls/ds/tools/forElogPost.txt"
    else:
        auth_path: str = f"/cds/home/opr/{opr_name}/forElogPost.txt"

    with open(auth_path, "r") as f:
        pw: str = f.readline()[:-1]

    return HTTPBasicAuth(username=opr_name, password=pw)

def getRunsWithTag(
        exp: str,
        tag: str,
        http_auth: Optional[HTTPBasicAuth] = None
) -> List[int]:
    """Return a list of runs tagged with a specific `tag`.

    Parameters
    ----------
    exp (str) Experiment name.
    tag (str) Tag to match against run tags.
    http_auth (HTTPBasicAuth) Authentication for eLog API.

    Returns
    -------
    tagged_runs (list[int]) List of runs with the specified tag. Empty if none
        were found or there was a communication error.
    """
    tag_url: str = f"{BASE_URL}/{exp}/ws/get_runs_with_tag?tag={tag}"
    http_auth: HTTPBasicAuth = http_auth or getElogBasicAuth(exp)
    resp: requests.models.Response = requests.get(tag_url, auth=http_auth)

    tagged_runs: list = []
    if resp.json()['success']:
        tagged_runs = resp.json()['value']

    return tagged_runs

def postElogMsg(
        exp: str,
        msg: str,
        *,
        tag: Optional[str] = "",
        title: Optional[str] = "",
        files: List[Union[str, Tuple[str, str]]] = []
) -> None:
    """Post a new message to the eLog. Adapted from `elog` package.

    Parameters
    ----------
    exp (str) Experiment name.
    msg (str) Body of the eLog post.
    tag (str) Optional. A tag to include for the post.
    title (str) Optional. A title for the eLog post.
    files (list) Optional. Either a list of paths (str) to files (figures) to
        include with the eLog post, OR, a list of 2-tuples of strings of the
        form (`path`, `description`).
    """
    post_files: List[Tuple[str, Tuple, Optional[str]]] = []
    for f in files:
        if isinstance(f, str):
            desc: str = os.path.basename(f)
            formatted_file: tuple = (
                "files",
                (desc, open(f, "rb")),
                mimetypes.guess_type(f)[0]
            )
        elif isinstance(f, tuple) or isinstance(f, list):
            formatted_file: tuple = (
                "files",
                (f[1], open(f[0], "rb")),
                mimetypes.guess_type(f[0])[0]
            )
        else:
            logger.debug(f"Can't parse file {f} for eLog attachment. Skipping.")
            continue
        post_files.append(formatted_file)

    post: Dict[str, str] = {}
    post['log_text'] = msg
    if tag:
        post['log_tags'] = tag
    if title:
        post['log_title'] = title

    http_auth: HTTPBasicAuth = getElogBasicAuth(exp)
    post_url: str = f"{BASE_URL}/{exp}/ws/new_elog_entry"

    params: Dict[str, Any] = {'url': post_url, 'data': post, 'auth': http_auth}
    if post_files:
        params.update({'files': post_files})

    resp: requests.models.Response = requests.post(**params)

    if resp.status_code >= 300:
        logger.debug(
            f"Error when posting to eLog: HTTP status code {resp.status_code}"
        )

    if not resp.json()['success']:
        logger.debug(f"Error when posting to eLog: {resp.json()['error_msg']}")

def postRunTable(exp: str, run: int, runtable_data: Dict[str, Union[int, float]]) -> None:
    ws_url: str = f"{BASE_URL}/run_control/{exp}/ws/add_run_params"
    http_auth: HTTPBasicAuth = getElogBasicAuth(exp)

    post_params: Dict[str, int] = {"run_num": run}

    resp: requests.models.Response = requests.post(
        ws_url,
        params=post_params,
        json=runtable_data,
        auth=http_auth
    )

    if resp.status_code >= 300:
        logger.debug(f"Error posting run table - status {resp.status_code}")

    if not resp.json()['success']:
        logger.debug(f"Error posting run table: {resp.json()['error_msg']}")
