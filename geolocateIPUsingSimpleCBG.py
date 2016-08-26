from cbg import cbg
from pprint import PrettyPrinter
import traceback
import json
import sys

if __name__=="__main__":
    try:
        pp=PrettyPrinter()
        cbgObj=cbg()


        measurementData=json.load(open(sys.argv[1],'r'))
        inputConstraints=cbgObj.getInputConstraints(measurementData)
        cityDict=cbgObj.getCities(inputConstraints)

        pp.pprint(cityDict)

    except:
        if len(sys.argv)==1:
            print('No input!')
        else:
            traceback.print_exc()

