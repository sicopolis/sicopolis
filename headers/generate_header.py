## Script to launch ensemble of simulations for one parameter


# one variable, more parameters

fldr = "/home/marius/Models/Sicopolis/headers"

file_in = "sico_specs_const_calv.h"


params = ["KONST_CALV"]

values = ["0.0d00","100.0d00","500d00","750d00","1000d00"]


for i in range(len(values)):

    #print(i)

    file_out = "sico_specs_const_calv_" + str(i+1) + ".h"

    with open(file_in, 'r') as fin, open(file_out, 'w') as fout:

        for line in fin:

            s = line.split()

            if line.startswith("#define") and s[1] in params:

                print(i)

                val = values[i]    

                line = "#define " + params[0] + " " + val + "\n"

            fout.write(line)
