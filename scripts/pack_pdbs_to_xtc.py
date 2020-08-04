
import biotite.structure as bs
import biotite.structure.io.pdb as bpdb
import biotite.structure.io.xtc as bxtc
import numpy as np
import argparse
import tarfile
import gzip
import pandas as pd
import os
import io
import time
import datetime

def load_pdbs_from_tar_gz(file_path,mode="r:gz"):
    """

    """
    if type(file_path) is not str:
        raise RuntimeError("Argument \'file_path\' is not of type str!")
    archive = tarfile.open(name=file_path,mode=mode)
    coords = []
    for pdb in archive.getmembers():
        if pdb.name == 'test_archive':
            continue
        pdb_handle = archive.extractfile(pdb)
        content = pdb_handle.read().decode()
        f = io.StringIO(content)
        pdbfile = bpdb.PDBFile()
        pdbfile.read(f)
        coords.append(pdbfile.get_coord()[0])
    archive.close()
    return np.array(coords)

def load_pdbs_from_dict(dir_path, sc=None):
    """

    """
    if type(dir_path) is not str:
        raise RuntimeError("Argument \'dir_path\' is not of type str!")
    files = os.listdir(dir_path)
    

    if sc is not None:
        print("using score")
        score = pd.read_csv(sc, header=1, delimiter= "\s+")         
        score = score[~pd.isnull(score['description'])]
        assert(len(score.keys()) > 0)
        assert('description' in score.keys())

#        for f in score['description']:
#            cond = f.split("/")[-1]+".pdb" in files
#            if not cond:
#                print(f)
#            assert(cond)            

        files = [os.path.join(dir_path,f.split("/")[-1]+".pdb") for f in score['description']]
    
    else:
        files = [os.path.join(dir_path,f) for f in files]

    coords = []
    i=0

    dT = 0

    N_files = len(files)


    #markers = np.linspace(0, 1, 100)
    p_marker_delta = 1.0/100

    marker = ''

    p_last = 0
    times = []

    print("first file ::" +str(files[0]))
    print("last  file ::" +str(files[-1]))
    
    for f in files:
        t1 = time.time()

        pdbfile = bpdb.PDBFile()
        pdbfile.read(f)

        try:
#            if coords is None:
#                coords = pdbfile.get_coord()[0]
#            else:
#                coords = np.vstack((coords, pdbfile.get_coord()[0]))

            coords.append(pdbfile.get_coord()[0])     
        except:

            print("")            
            print("")
            print(" tried loading file ::"  + str(f))   
            print("")            
            exit()
        
        t2 = time.time()

        times.append(t2-t1)
        
        p = float(i+1)/float(N_files)
        
        if p - p_last >= p_marker_delta:
            p_last = p
            marker = marker + "="    

        T_eta = np.mean(times)*(N_files-(i+1))
        T_eta_err = np.std(times)

        print(
            "loading file ["  +str("") + " | " + str(i) + " / " + str(N_files) + "]  " +
             "|" + str(marker) + ">   " + str(p*100)+ "% complete" +
             " ETA :: " + str(datetime.timedelta(seconds=T_eta)) + " +/- " + str(T_eta_err) + "s ", 
            end='\r'
        )        

        i = i+1
    print("")
    print("")
    print(" ... finished loading pdbs ... ")

    return np.array(coords, dtype=np.float64)

def load_pdbs(files, sc=None):
    """

    """
    if type(files) is not str:
        raise RuntimeError("Argument \'files\' is not of type str!")
    endings = files.split(".")
    if "tar" in endings:
        if "gz" in endings:
            stack = load_pdbs_from_tar_gz(files,"r:gz")
        else:
            stack = load_pdbs_from_tar_gz(files,"r")
    elif os.path.isdir(files):
        stack = load_pdbs_from_dict(files, sc)
    else:
        raise RuntimeError("\'" + files + "\' is neither a tar archive nor a directory. Check your input.")
    return stack

def write_to_xtc(array_stack,outfile):
    """
    """
    xtcfile = bxtc.XTCFile()
    xtcfile.set_coord(array_stack)
    xtcfile.write(outfile)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("convert pdb files in directory or tar archive into xtc file.")
    parser.add_argument('-i',action='append')
    parser.add_argument('-o',action='append')
    parser.add_argument('-sc',action='append')
    args = parser.parse_args()
    if args.sc is None:
        stack = load_pdbs(files = args.i[0])
    else:
        stack = load_pdbs(files = args.i[0], sc = args.sc[0])
    write_to_xtc(array_stack = stack, outfile = args.o[0])

