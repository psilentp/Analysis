import db_access as dba
import pylab as plb

db = dba.get_db()

def cefig(i,cenum):
    f = plb.figure(i,figsize = (3,6))
    db[cenum].long_ifam.plot()
    f.axes[0].set_ybound(-50,190)
    f.axes[0].set_xbound(0,0.04)
    plb.show()

if __name__ == '__main__':
    AVBID = 1374
    AVEID = 1325
    AVAID = 1331
    AIYID = 1252
    
    for i,cennum in enumerate([AVAID,AVBID,AVEID,AIYID]):
        cefig(i,cennum)

