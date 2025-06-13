#!/usr/bin/env python
import os
import numpy as np
import argparse
import logging
from scipy import sparse
from pathlib import Path
from PIL import Image

def get_detinfo(fh5, detname):
    imgShape = getattr(getattr(fh5.UserDataCfg, detname), 'imgShape').read()
    ix = getattr(getattr(fh5.UserDataCfg, detname), 'ix').read()
    iy = getattr(getattr(fh5.UserDataCfg, detname), 'iy').read()
    pixel_array = np.ones_like(ix)

    mask = getattr(getattr(fh5.UserDataCfg, detname), 'mask').read()
    if (np.array(mask.shape)-np.array(imgShape)).sum()!=0:
        piximg = np.asarray(
            sparse.coo_matrix(
                (pixel_array.flatten(), (ix.flatten(), iy.flatten())), shape=imgShape
            ).todense()
        )
    else:
        piximg = np.ones_like(mask)
        ret_dict={'ix':ix,
              'iy':iy,
              'imgShape':imgShape,
              'piximg':piximg}
    return ret_dict

# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define Args
parser = argparse.ArgumentParser()
parser.add_argument(
    "--run", help="run", type=str, default=os.environ.get("RUN_NUM", "")
)
parser.add_argument(
    "--experiment",
    help="experiment name",
    type=str,
    default=os.environ.get("EXPERIMENT", ""),
)
parser.add_argument(
    "--sum",
    help="save sum images as tiff",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--sumavg",
    help="save average images as tiff",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--events",
    help="save tiff files for single events",
    type=str,
    default=''
)
parser.add_argument(
    "--directory",
    help="directory to read files from (def <exp>/hdf5/smalldata)",
    default=None,
)
args = parser.parse_args()
logger.debug("Args to be used for pedestal plots: {0}".format(args))

expname = args.experiment
run = int(args.run)
S3DF_BASE = Path("/sdf/data/lcls/ds/")
hutch = expname[:3]
dirname = (
    f"{S3DF_BASE}/{hutch}/{expname}/hdf5/smalldata"
)
if args.directory: dirname = args.directory
fname = "%s/%s_Run%04d.h5" % (dirname, expname, run)
import tables
if os.path.isfile(fname):
    fh5 = tables.open_file(fname).root
else:
    print(fname,' done not exist, exit')
    sys.exit()

tiffdirname = S3DF_BASE / f"{expname[:3]}/{expname}/scratch/run{int(run)}"
if not os.path.isdir(tiffdirname):
    os.makedirs(tiffdirname)

alldetinfo={}
    
if args.sum or args.sumavg:
    nEvts = fh5.event_time.shape[0]
    #now I need to get keys....
    sumnames = [dmg for dmg in dir(getattr(fh5, 'Sums')) if dmg[0]!='_']

    for sumname in sumnames:
        data = getattr(fh5, f"Sums/{sumname}").read().squeeze()
        detname = sumname.split('_calib')[0]
        if detname in alldetinfo:
            detinfo = alldetinfo[detname]
        else:
            detinfo = get_detinfo(fh5, detname)
            alldetinfo[detname]=detinfo
        if (detinfo['imgShape']-data.shape).sum()!=0:
            img = np.asarray(
                sparse.coo_matrix(
                    (data.flatten(),
                     (detinfo['ix'].flatten(), detinfo['iy'].flatten())),
                    shape=detinfo['imgShape']
                ).todense())
            data = img/detinfo['piximg']
        data = data/nEvts

        im = Image.fromarray(data)
        tiff_file = f"{tiffdirname}/Run_{int(run)}_{sumname}.tiff"
        print(tiff_file)
        im.save(tiff_file)

            
if args.events != "":
    print('Looking at specified detectors only: ',args.events)
    detnames = args.events.split(',')

    for detname in detnames:
        if detname in alldetinfo:
            detinfo = alldetinfo[detname]
        else:
            detinfo = get_detinfo(fh5, detname)
            alldetinfo[detname]=detinfo

        alldata = [d for d in dir(getattr(fh5, f"{detname}")) if d[0]!='_']
        for d in alldata:
            shp = getattr(getattr(fh5, f"{detname}"),d).shape
            if len(shp)>2:
                for evt in range(shp[0]):
                    data = getattr(getattr(fh5, f"{detname}"),d).read(evt,evt+1).squeeze()
                    if ((detinfo['imgShape']-data.shape).sum())!=0:
                        img = np.asarray(
                            sparse.coo_matrix(
                                (data.flatten(),
                                 (detinfo['ix'].flatten(), detinfo['iy'].flatten())),
                                shape=detinfo['imgShape']
                            ).todense())
                        data = img/detinfo['piximg']
                        
                    im = Image.fromarray(data)
                    tiff_file = f"{tiffdirname}/Run_{int(run)}_evt_{evt+1}_{detname}.tiff"
                    im.save(tiff_file)

    #fh5.close()
