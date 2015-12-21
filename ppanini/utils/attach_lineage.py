import os
import sys
import pdb
import re
'''Tmp file to parse results'''
def read_lineage(filename):
	lineage = {}
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]', '', i) for i in line.split('|')]
			try:
				g_ind = [i for i in range(len(split_line)) if 'g__' in split_line[i]][0]
			except:
				# pdb.set_trace()
				raise Exception('no g__ level found in hierarchy'+line)
			
			lineage['.'.join(split_line[g_ind:g_ind+2])] = '.'.join(split_line)
	return lineage

def sim_lineage(filename):
	
	s_g, g_f, f_o, o_c, c_p, p_k ={}, {}, {}, {}, {}, {}

	lineage = {'s_g':s_g, \
			   'g_f':g_f, \
			   'f_o':f_o, \
			   'o_c':o_c, \
			   'c_p': c_p, \
			   'p_k':p_k}

	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]', '', i) for i in line.split('|')]

			s_ind = [i for i in range(len(split_line)) if 's__' in split_line[i]]
			if s_ind:
				s_ind = s_ind[0]
			
			g_ind = [i for i in range(len(split_line)) if 'g__' in split_line[i]]
			if g_ind:
				g_ind = g_ind[0]

			f_ind = [i for i in range(len(split_line)) if 'f__' in split_line[i]]
			if f_ind:
				f_ind = f_ind[0]
			o_ind = [i for i in range(len(split_line)) if 'o__' in split_line[i]]
			if o_ind:
				o_ind = o_ind[0]
			c_ind = [i for i in range(len(split_line)) if 'c__' in split_line[i]]
			if c_ind:
				c_ind = c_ind[0]
			p_ind = [i for i in range(len(split_line)) if 'p__' in split_line[i]]
			if p_ind:
				p_ind = p_ind[0]
			k_ind = [i for i in range(len(split_line)) if 'k__' in split_line[i]]
			if k_ind:
				k_ind = k_ind[0]

			if s_ind:
				if not split_line[s_ind] in s_g:
					s_g[split_line[s_ind]] = split_line[g_ind]
			if g_ind:
				if not split_line[g_ind] in g_f:
					g_f[split_line[g_ind]] = split_line[f_ind]
			if f_ind:
				if not split_line[f_ind] in f_o:
					f_o[split_line[f_ind]] = split_line[o_ind]
			if o_ind:
				if not split_line[o_ind] in o_c:
					o_c[split_line[o_ind]] = split_line[c_ind]
			if c_ind:
				
				if not split_line[c_ind] in c_p:
					c_p[split_line[c_ind]] = split_line[p_ind]
			if p_ind:
				
				if not split_line[p_ind] in p_k:
					p_k[split_line[p_ind]] = split_line[k_ind]
	return lineage

if __name__ == "__main__":
	
	help = ['-h', '--help', '--h', '-help']
	if  sys.argv[1] in help:
		sys.exit('python '+sys.argv[0]+' <filename> <taxonomy> > <output_file>')

	lineage = read_lineage(sys.argv[2])
	simulated_lineage = sim_lineage(sys.argv[2])
	filename = sys.argv[1]
	sim_lineage_bool = False
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')]
			g_s = [re.sub('[\r\t\n]','', i) for i in split_line[-1].split('.')]
			g_ind = [i for i in range(len(g_s)) if 'g__' in g_s[i]][0]

			genome_name = '.'.join(g_s[g_ind:g_ind+2])
			if genome_name in lineage:
				print split_line[0]+'\t'+lineage[genome_name]
			else:
				if g_s[g_ind] in simulated_lineage['g_f']:
					genus_i = g_s[g_ind]
					family_i = simulated_lineage['g_f'][genus_i]
					order_i = simulated_lineage['f_o'][family_i]
					class_i = simulated_lineage['o_c'][order_i]
					phylum_i = simulated_lineage['c_p'][class_i]
					kingdom_i = simulated_lineage['p_k'][phylum_i]
					
					lineage_i = '.'.join([kingdom_i, phylum_i, class_i, order_i, family_i]+g_s)
					print split_line[0]+'\t'+lineage_i
				else:
					print 'taxonomy NOT COMPLETE: '+split_line[-1]+' NOT FOUND'	

