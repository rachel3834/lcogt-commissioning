'''
Tool to copy gain history data from one database to another, in particular to support hte migration from local sqlite to online postgress.
'''
import argparse
import logging

from lcocommissioning.common.noisegaindb_orm import noisegaindb, NoiseGainMeasurement


def parseCommandLine():
    parser = argparse.ArgumentParser(
    description='Copy noisegain database from A to B')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                    help='Set the debug level')

    parser.add_argument('--inputurl', type=str,  default='sqlite:///noisegain.sqlite', help="input database")
    parser.add_argument('--outputurl' ,type=str , default='sqlite:///noisegain.sqlite', help="input database")

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args

if __name__ == '__main__':
    args = parseCommandLine()
    print (f"Copy from {args.inputurl} -> {args.outputurl}")

    input = noisegaindb(args.inputurl)
    output = noisegaindb (args.outputurl)


    q = input.session.query (NoiseGainMeasurement)
    print ("Found {} records to copy.".format(q.count()))
    newdata = [NoiseGainMeasurement(e) for e in q.all()]

    print ("Now doing bulk insert")
    output.session.bulk_save_objects (newdata)

    output.session.commit()

    input.close()
    output.close()
