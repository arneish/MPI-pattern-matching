#!/usr/bin/python3

#########################################################################
# Generate the input file containing text and pattern in required       #
# format and save to file named testcase_<n>_<num_patterns>             #
#                                                                       #
# Parameters:                                                           #
#   n               :length of text                                     #
#   num_patterns    :number of patterns to be searched                  #
#   min_p           :minimum period length                              #
#   min_m           :minimum period length                              #
# Format of output file:                                                #
#   -----------------------------------------------------------------   #
#   | n num_patterns                                                    #
#   | text                                                              #
#   | m[0] m[1] m[2] ... m[num_patterns-1]                              #
#   | p[0] p[1] p[2] ... p[num_patterns-1]                              #
#   | pattern[0]                                                        #
#   | pattern[1]                                                        #
#   | ...                                                               #
#   | pattern[num_pattern-1]                                            #
#   -----------------------------------------------------------------   #
#                                                                       #
#########################################################################

import random
import string

def randomString(stringLength):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

# adjust these parameters to generate testcases
n            = 10000    # length of text
num_patterns = 10       # number of patterns to be searched
min_p        = 2        # minimum period length
min_m        = 5        # minimum pattern length


m_set = []
p_set = []
pattern_set = []

filename = 'testcase_' + str(n) + '_' + str(num_patterns)     #output filename
file = open(filename, 'w')

# write size of length of text and number of patterns in first line of file
file.write(str(n) + ' ' +str(num_patterns) + '\n')

# generate and write text
file.write(randomString(n) + '\n')

# generate and write length of patterns
for i in range(num_patterns):
	m = random.randint(min_m, n)
	m_set.append(m)
	file.write(str(m) + ' ')
file.write('\n')

# generate and write period length of patterns
for i in range(num_patterns):
	p = random.randint(min_p, m_set[i]//2)
	p_set.append(p)
	file.write(str(p) + ' ')
file.write('\n')

# generate and write patterns
for i in range(num_patterns):
	period_str = randomString(p_set[i])
	pattern = period_str * (m_set[i] // p_set[i])
	pattern = pattern + period_str[0 : (m_set[i]-len(pattern))]
	file.write(pattern + '\n')

file.close()
