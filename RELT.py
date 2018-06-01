"""
Usage: estimate heritability by regression of sample phenotypic covariance matrix onto genetic relatedness

This script performs relatedness thresholded (RELT) heritability estimation. It regresses elements of the sample
phenotypic covariance matrix onto corresponding elements of a relatedness matrix for those pairs with relatedness
below a user-set threshold (defaul 0.05).

The script calculates standard errors of genetic variance and heritability estimates by a procedure that takes into
account dependence between pairs. This can be computationally demanding for large sample sizes.

It takes a binary relatedness matrix as input. The matrix is formatted the same way as produced by
GCTA with the --make-grm-bin option set. The matrix is in lower-triangular order with 32 bit floating point numbers.

Associated with the relatedness matrix is an plain text id file. The first column of the id file gives the ids of the
individuals in the order that they appear in the relatedness matrix. If IDs are specified uniquely by the first column
of the GCTA grm.id file, then the GCTA grm.id file can be used.

The trait file is a plain text file with columns: FID, IID, trait1, trait2, etc.

If covariates are supplied, the trait will first be adjusted for covariates before heritability is estimated.

The script outputs a file outprefix.herit that records variance component estimates and standard errors.
"""
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ids', type=str, help='Location of file with ids in order of relatedness matrices')
parser.add_argument('R', type=str, help='Location of file giving the relatedness matrix')
parser.add_argument('trait', type=str, help='Location of file with trait file')
parser.add_argument('outprefix', type=str, help='Output prefix')
parser.add_argument('--rel_thresh',type=float,help='Maximum relatedness of pairs used to estimate heritability',default = 0.05)
parser.add_argument('--covariates', type=str, help='File with covariate matrix', default=None)
parser.add_argument('--trait_index', type=int, help='Index of trait to analyse in trait file', default=1)
args = parser.parse_args()


def rd_calc(rdp):
    # compute relatedness disequilibrium
    print('Computing regression estimate')
    xtx = np.zeros((2, 2), dtype=np.float64)
    for i in xrange(0, 2):
        for j in xrange(0, i + 1):
            xtx[i, j] = np.sum(rdp[i] * rdp[j])
            xtx[j, i] = xtx[i, j]
    xtx_inv = np.linalg.inv(xtx)
    Rd = np.zeros(rdp[0].shape, dtype=np.float32)
    for i in xrange(0, 2):
        print(str(xtx_inv[1, i]))
        Rd += xtx_inv[1, i] * rdp[i]

    return Rd


def read_trait(trait_file):
    trait_f = open(trait_file, 'r')
    trait_line = trait_f.readline()
    n_cols = len(trait_line.split(' '))
    trait_f.close()
    coltypes = ['S10', 'S10']
    for i in xrange(0, n_cols - 2):
        coltypes.append('f8')
    trait = np.genfromtxt(trait_file, dtype=coltypes, delimiter=' ', missing_values='NA')
    return trait


### Read relatedness matrix IDs ###
ids = np.loadtxt(args.ids, dtype='S10')
if len(ids.shape) > 1:
    ids = ids[:, 0]
N = len(ids)
id_dict = {}
for i in xrange(0, N):
    id_dict[ids[i]] = i
print('IDs of ' + str(N) + ' individuals read for relatedness matrices\n')

### Read trait ###
# Find number of columns in trait file

trait = read_trait(args.trait)
trait_ids = np.array([x[0] for x in trait], dtype='S10')
trait = np.array([x[args.trait_index + 1] for x in trait], dtype=np.float64)
not_na = np.logical_not(np.isnan(trait))
trait = trait[not_na]
trait_ids = trait_ids[not_na]
n_trait = np.sum(not_na)

print('Non-missing phenotypes of ' + str(n_trait) + ' individuals included')

trait_dict = {}
for i in xrange(0, n_trait):
    trait_dict[trait_ids[i]] = i

############# Adjust for covariates #############
if args.covariates is None:
    trait = trait - np.mean(trait)
else:
    # Read covariates
    covars = read_trait(args.covariates)
    ncovar = len(covars[0]) - 2
    covar_n = covars.shape[0]
    print(str(ncovar) + ' covariates included')
    # Match covariate ids to trait ids
    covar_ids = np.array([x[0] for x in covars], dtype='S10')
    covar_id_set = set(covar_ids)
    # Remove trait observations that do not have covariate observations
    trait_ids_in_covar = np.array([x in covar_id_set for x in trait_ids])
    trait_ids = trait_ids[trait_ids_in_covar]
    trait = trait[trait_ids_in_covar]
    n_trait = trait.shape[0]
    trait_dict = {}
    for i in xrange(0, n_trait):
        trait_dict[trait_ids[i]] = i
    # Find covariate ids that are in trait ids
    covar_id_in_trait = np.array([x in trait_dict for x in covar_ids])
    print(str(n_trait) + ' individuals in common between trait and covariate observations')
    # Make covariate matrix matching trait
    covar_matrix = np.ones((n_trait, ncovar + 1), dtype=np.float64)
    for i in xrange(0, covars.shape[0]):
        if covar_id_in_trait[i]:
            for j in xrange(0, ncovar):
                covar_matrix[trait_dict[covar_ids[i]], 1 + j] = covars[i][j + 2]
    # Adjust for covariates
    cxtx = np.dot(covar_matrix.T, covar_matrix)
    cxtx = np.linalg.inv(cxtx)
    trait_pred = np.dot(covar_matrix.dot(cxtx), np.dot(covar_matrix.T, trait))
    trait = trait - trait_pred

################# match up trait IDs and relatedness IDs ########

# Common between IDS and trait
ids_in_trait = np.zeros((N), dtype=bool)
for i in xrange(0, N):
    if ids[i] in trait_dict:
        ids_in_trait[i] = True

n = np.sum(ids_in_trait)
if n == 0:
    raise (ValueError('No individuals with both relatedness and trait information'))
print(str(n) + ' individuals with both relatedness and trait information')

trait_indices = np.zeros((n), dtype=int)
count = 0
for i in xrange(0, N):
    if ids_in_trait[i]:
        trait_indices[count] = trait_dict[ids[i]]
        count += 1

trait = trait[trait_indices]
trait_ids = trait_ids[trait_indices]


########### Read in relatedness matrices ##########
def grm_read(grm_file, ids_in_trait):
    n = np.sum(ids_in_trait)
    N = len(ids_in_trait)
    print('Reading ' + grm_file)
    R = np.zeros((N, N), dtype=np.float32)
    R[np.tril_indices(N)] = np.fromfile(grm_file, dtype=np.float32)
    R = R + R.T
    R[np.diag_indices(N)] = np.diag(R) / 2.0
    R = R[np.ix_(ids_in_trait, ids_in_trait)]
    return (R)


R_file = args.R
R = grm_read(R_file, ids_in_trait)

############ Estimate phenotypic covariance ######
diag_R = np.zeros((n))
diag_R[:] = np.diag(R)
# set diagonal to zero
np.fill_diagonal(R, 0)

# constant matrix for regression
S = np.ones(R.shape, dtype=np.float32)
np.fill_diagonal(S, 0)

# Compute regression matrix
R_reg = rd_calc([S, R])

v = np.dot(trait.T, R_reg.dot(trait))

del R_reg

# print('v_g: '+str(np.round(v_g,3))+' ('+str(np.round(np.sqrt(v_g_var),3))+')')
print('v: ' + str(np.round(v, 3)))

v_y = np.sum(np.square(trait)) / np.float(n)
print('Phenotypic variance estimate: ' + str(np.round(v_y, 3)))

R_mean = np.sum(R) / (n ** 2 - n)
sigma = v * (R - R_mean) + (v_y - v) * np.diag(np.ones((n), dtype=np.float32))
sigma[np.diag_indices(n)] = (v_y - v) + v * diag_R

################ Estimate heritability (thresholded) ##########
print('Estimating heritability')
close_rel = R > args.rel_thresh
n_close_rel = np.sum(close_rel)
print(str(n_close_rel) + ' pairs with relatedness >'+str(args.rel_thresh)+' excluded')
R[close_rel] = 0
S[close_rel] = 0
del close_rel
R_reg = rd_calc([S, R])
del R
v_g = np.dot(trait.T, R_reg.dot(trait))

print('v_g: ' + str(np.round(v_g, 3)))

print('Estimating standard error')

sigma_R = sigma.dot(R_reg)
v_g_var = 2 * np.sum(np.square(sigma_R))

print('v_g: ' + str(np.round(v_g, 3)) + ' (' + str(np.round(np.sqrt(v_g_var), 3)) + ')')

print('Estimating variance of phenotypic variance estimate')
var_v_y = 2 * np.sum(np.square(sigma)) / np.float(n ** 2)

print('Estimating covariance between genetic variance and phenotypic variance estimate')

cov_v_g_v_y = np.sum(sigma * sigma_R)

cov_v_g_v_y = 2 * cov_v_g_v_y / np.float(n)

h2 = v_g / v_y

h2_var = ((v_g ** 2) / (v_y ** 2)) * (v_g_var / (v_g ** 2) - 2 * cov_v_g_v_y / (v_g * v_y) + var_v_y / (v_y ** 2))

print('v_g/v_y: ' + str(np.round(h2, 3)) + ' (' + str(np.round(np.sqrt(h2_var), 4)) + ')')

herit_out = open(args.outprefix + '.herit', 'w')
herit_out.write('\tEstimate\tS.E.\n')
herit_out.write('v\t' + str(v) + '\tNA\n')
herit_out.write('v_g\t' + str(v_g) + '\t' + str(np.sqrt(v_g_var)) + '\n')
herit_out.write('v_y\t' + str(v_y) + '\t' + str(np.sqrt(var_v_y)) + '\n')
herit_out.write('h2\t' + str(h2) + '\t' + str(np.sqrt(h2_var)) + '\n')
herit_out.close()
