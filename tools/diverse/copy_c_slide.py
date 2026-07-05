import sys

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

'''
Copy an entry for C_SLIDE_DIMLESS in a SICOPOLIS run-specs header 
to a different (existing) header.

Execution of the script from the main SICOPOLIS directory:
  python3.11 ./tools/diverse/copy_c_slide.py <run_name_1> <run_name_2>
(A newer version of Python will also do.)

Run-specs headers must be in the standard directory './headers';
otherwise, the script won't work!

Author: Ralf Greve
Last update: 2026-07-05
'''

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------- Run-specs headers --------

if len(sys.argv) < 3:
    print('Error: Missing <run_name> arguments.')
    print('Usage: python copy_c_slide.py <run_name_1> <run_name_2>.')
    sys.exit(1)  # Stop the script immediately

run1 = sys.argv[1]
run2 = sys.argv[2]

header_file_1 = f'sico_specs_{run1}.h'
header_file_2 = f'sico_specs_{run2}.h'

#-------- Read C_SLIDE_DIMLESS from first header --------

with open(f'./headers/{header_file_1}', 'r', encoding='ascii') as header1:
    content1 = header1.readlines()

line_num_1 = None
for num, line in enumerate(content1):
    if '#define C_SLIDE_DIMLESS' in line:
        line_num_1 = num
        break
if line_num_1 is None:
    print(f'Error: \'C_SLIDE_DIMLESS\' was not found in the header1 file!')
    exit()

c1 = content1[line_num_1]

#-------- Write C_SLIDE_DIMLESS to second header,
#                       replacing previous entry --------

with open(f'./headers/{header_file_2}', 'r', encoding='ascii') as header2:
    content2 = header2.readlines()

line_num_2 = None
for num, line in enumerate(content2):
    if '#define C_SLIDE_DIMLESS' in line:
        line_num_2 = num
        break
if line_num_2 is None:
    print(f'Error: \'C_SLIDE_DIMLESS\' was not found in the header2 file!')
    exit()

replacement = c1 if c1.endswith('\n') else c1 + '\n'

content2[line_num_2] = replacement

with open(f'./headers/{header_file_2}', 'w', encoding='ascii') as header2:
    header2.writelines(content2)

#-------- End of script --------

print('Done.')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
