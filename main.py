__author__ = 'psilentp'
from read_heka import BundleHeader,PGFFile,PULFile, getleafs
def main():
    import neo
    f = open('./test_data/CEN184/THL_2012-03-21_18-40-42_000.dat')
    #create a bundle header object with the file
    head = BundleHeader(f)
    #load the file
    head.load(f)
    #get the .pgf and .pul items in the file
    for bi in head.oBundleItems:
        if str(bi.oExtension)[0:4] == '.pgf':
            pgf = PGFFile(f,bi)
        if str(bi.oExtension)[0:4] == '.pul':
            pul = PULFile(f,bi)
    return pul

if __name__ == '__main__':
    main()