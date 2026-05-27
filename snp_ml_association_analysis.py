import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy as sp
from bioinfokit import analys, visuz

# Read the data
temp_genotypes = pd.read_excel('genotypes.xls', index_col=0, skiprows=1)
temp_phenotypes = pd.read_excel('phenotypes.xls', index_col=0)

# define Phenotype as index
temp_phenotypes = temp_phenotypes.set_index('Phenotype')

# search for my phenotype
phenotype_name = 'Tumor necrosis factor (TNF)-alpha level in lung after ' \
                 'aerosolized lipopolysaccharide (LPS) exposure [pg/ml]'
my_phenotype = temp_phenotypes.loc[phenotype_name]

# clean my phenotype data from null values
phenotype = my_phenotype.dropna()

# clean genotype data from unknown alleles
genotypes = temp_genotypes.drop(['BXD101', 'BXD102', 'BXD103'], axis=1)

# save the information about each snp
snp_details = genotypes[['Locus', 'Chr_Build37', 'Build37_position']]

# define locus as index
genotypes = genotypes.set_index('Locus')

# clean the table
genotypes = genotypes.drop(['Chr_Build37', 'Build37_position'], axis=1)

# convert all letters to uppercase
genotypes = genotypes.apply(lambda allele: allele.str.upper(), axis=1)


# Statistical F test
def f_test(r_square, n, k):
    dfn = n - 2
    dfd = k - 1
    f_statistics = r_square / ((1 - r_square) / (n - 2))
    # print("F test: ", f_statistics)
    p_val = 1 - sp.stats.f.cdf(f_statistics, dfd, dfn)
    return p_val


# Linear regression function
def linear_reg(x, y):
    k = 2
    # print("k: ", k)
    n = x.shape[0]
    # print("n: ", n)
    mean_x = x.mean()
    # print("mean x: ", mean_x)
    mean_y = y.mean()
    # print("mean y: ", mean_y)
    var_x = x.var()
    # print("var_x: ", var_x)
    cov_x_y = (sum(x * y) - (n * mean_x * mean_y)) / (n - 1)
    # print("cov (x,y): ", cov_x_y)
    # calculate b0, b1
    b1 = cov_x_y / var_x
    b0 = mean_y - (b1 * mean_x)
    # calculate p-value
    y_hat = b0 + b1 * x
    # print("y_hat: ", y_hat)
    sse = sum(np.square(y - y_hat))
    # print("SSE: ", sse)
    ssr = sum((np.square(y_hat - mean_y)))
    # print("SSR: ", ssr)
    sst = sse + ssr
    # print("SST: ", sst)
    r_square = ssr / sst
    # print("r square: ", r_square)
    p_val = f_test(r_square, n, k)
    return b0, b1, p_val


# Loop over all SNPs and calculate the -log(p-val)
p_val_dict = dict()
for snp, snp_info in zip(genotypes.iterrows(), snp_details.iterrows()):
    snp_without_H = snp[snp != 'H']
    # keep only columns that exists in both tables
    phenotype_align, snp_align = phenotype.align(snp_without_H, join="inner", axis=0)
    snp_align = snp_align.apply(
        lambda strain: 0 if strain == 'B' else 2)  # change to numeric values
    p_val = linear_reg(snp_align.to_numpy(), phenotype_align.to_numpy())[2]
    p_val_dict[snp_info[1]['Locus']] = [p_val, snp_info[1]['Chr_Build37']]

# create p_value df for output (contains position, chromosome and -log(p-value) result for each SNP)
p_values_output = pd.DataFrame.from_dict(p_val_dict, orient='index')
p_values_output.rename(columns={0: 'p-value', 1: 'Chr'}, inplace=True)

# creating Manhattan Plot
p_val_to_plot = p_values_output.copy()
visuz.marker.mhat(df=p_values_output, chr='Chr', pv='p-value')

# Correct for multiple compression
corrected_alpha = 0.05 / (p_values_output.shape[0])
p_values_output[p_values_output['p-value'] <= corrected_alpha]

(p_values_output.drop(columns=['p-value']).rename(columns={'tpval': '-log10(p-value)'})).to_excel(
    '-log_p_values_output.xlsx')
