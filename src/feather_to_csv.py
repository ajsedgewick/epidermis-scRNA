#!/usr/bin/env python3

import sys
import feather

def main(args):

    try:
        df = feather.read_dataframe(args[0])
        df.set_index(df.iloc[:,0], inplace=True)
        df.index.name=None
        df.drop(df.columns[0],axis=1,inplace=True)
        df.to_csv(args[1], float_format = "%.1f")
    except OSError as err:
        print("OS error: {0}".format(err))

        
if __name__ == "__main__":
    main(sys.argv[1:])
