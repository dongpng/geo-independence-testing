#!/usr/bin/env python
import optparse
import sys

import numpy as np

import HSIC

def main():
    # Check options
    p = optparse.OptionParser()
    p.add_option('--datatype', '-d', help="The type of data: bin (binary), cat (categorical), freq (frequency)")
    p.add_option('--locationsfile', '-l', help="The file with the locations")
    p.add_option('--featurefile', '-f', help="The file with the linguistic data")
    p.add_option('--num_samples', '-n', help="Number of bootstrap samples for calculating p-values", default=500)

    options, arguments = p.parse_args()
    if not options.locationsfile:   
        p.error('Inputfile with locations is not given. Each line is an observation (coordinates are seperated by a space).')
    if not options.featurefile:   
        p.error('Inputfile with linguistic data is not given. Each line is an observation.')
    if not options.datatype:   
        p.error('Datatype (bin/cat/freq) is not given')

    # Read in data
    locs = []
    with open(options.locationsfile, 'r') as input_file:
        for line in input_file:
            values = line.split(" ")
            locs.append([float(values[0]), float(values[1])])

    locs = np.array(locs)

    data = []
    with open(options.featurefile, 'r') as input_file:
        for line in input_file:
            data.append([float(line)])

    data = np.array(data)   
   
    # Check if they are the same size
    if len(data) != len (locs):
        print "Length of datafiles are not equal!"
        sys.exit(0)


    if options.datatype == "freq":
        print("HSIC: %.5f (p-value: %.3f)" % HSIC.HSIC_pval(locs, data, 
                            kernelX="Gaussian", kernelY="Gaussian", N_samp=options.num_samples))

    elif options.datatype == "bin":
        print("HSIC: %.5f (p-value: %.3f)" % HSIC.HSIC_pval(locs, data,
                            kernelX="Gaussian", kernelY="Delta", N_samp=options.num_samples, p_method="boots")) 

    elif options.datatype == "cat":
        print("HSIC: %.5f (p-value: %.3f)" % HSIC.HSIC_pval(locs, data, 
                            kernelX="Gaussian", kernelY="Delta", N_samp=options.num_samples, p_method="boots")) 




if __name__ == '__main__':
    main()

   