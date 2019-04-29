import hddm
import matplotlib.pyplot as plt

# Load data from csv file into a NumPy structured array
data = hddm.load_csv('performance_all_LSB.csv')
print(data.head())

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data.groupby('subj_idx'):
    subj_data.rt.hist(bins=20, histtype='step', ax=ax)
plt.show()

# Instantiate model object passing it our data (no need to call flip_errors() before passing it).
# This will tailor an individual hierarchical DDM around your dataset.
m = hddm.HDDM(data)
# find a good starting point which helps with the convergence.
m.find_starting_values()
# start drawing 7000 samples and discarding 5000 as burn-in
m.sample(2000, burn=20)

stats = m.gen_stats()
stats[stats.index.isin(['a', 'a_std', 'a_subj.0', 'a_subj.1'])]
m.plot_posteriors(['a', 't', 'v', 'a_std'])
plt.show()

m.plot_posterior_predictive(figsize=(14, 10))
plt.show()

# Create a HDDM model multi object
#model = hddm.HDDM(data, depends_on={'v':'difficulty'})

# Create model and start MCMC sampling
#model.sample(2000, burn=20)

# Print fitted parameters and other model statistics
#model.print_stats()

# Plot posterior distributions and theoretical RT distributions
#model.plot_posteriors()
#plt.show()
#model.plot_posterior_predictive()
#plt.show()
