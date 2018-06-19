#!/usr/bin/env python3

import sys
import feather
import pandas as pd

def main(args):

    try:
        df = pd.read_csv(args[0], index_col=0, header=0)
        df.insert(loc=0, column="rowname", value=df.index)
        feather.write_dataframe(df, args[1])

    except OSError as err:
        print("OS error: {0}".format(err))

        
if __name__ == "__main__":
    main(sys.argv[1:])
