import mpc
import pylab as lab
import sys

# make sure root is startet in Batch mode
sys.argv.append('-b')
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.StartGuiThread = False 

gamma = lab.zeros(200)
for iE in xrange(200):
    gamma[iE] = 10. ** (6. + (iE * 8. / 199.))


def get_rate_h(interaction, id, gamma, count=1000):
    h = {}
    c = mpc.Candidate()
    c.current.setId(id)
    c.current.setLorentzFactor(gamma)
    c.setCurrentStep(0.0 * mpc.Mpc)

    for i in range(count):
        c.setNextStep(10000000 * mpc.Mpc)
        c.clearInteractionStates()
    
        interaction.process(c)
        state = mpc.InteractionState()
        c.getInteractionState("mpc::PhotoDisintegration", state)
        if not h.has_key(state.channel):
            h[state.channel] = ROOT.TH1F('', '', 100, 0, 1000)
        h[state.channel].Fill(state.distance / mpc.Mpc)
    return h

def get_rates_per_gamma(interaction, id, gamma):
    r = {}
    
    h = get_rate_h(interaction, id, gamma, 10000)
    for key in h.keys():
        if h[key].GetEffectiveEntries() == 0:
            r[key] = 0
            continue
        f = ROOT.TF1('f1', 'expo')
        h[key].Fit(f, "q")
        r[key] = -f.GetParameter(1)
    return r

def extract_ids(table):    
    ids = []
    for i in xrange(len(table)):
        id = int(table[i][0])
        if not id in ids:
            ids.append(id)
    return ids

def get_rates_per_channel(interaction, id):
    rates = {}
    for iE in xrange(200):
        r = get_rates_per_gamma(interaction, id, gamma[iE])
        for channel in r.keys():
            if not rates.has_key(channel):
                rates[channel] = lab.zeros(200)
            rates[channel][iE] = r[channel]
    return rates

def get_data_rates_per_channel(table, id):
    rates = {}
    for i in xrange(len(table)):
        if int(table[i][0]) != id:
            continue
        rates[table[i][1]] = table[i][2:]
    return rates

def plot_channel(rates_simulated, rates_data, id, channel):
    if channel in rates_data:
        lab.plot(gamma, rates_data[channel], 'r', label="data " + str(channel), linewidth=2, markeredgewidth=2)
    if channel in rates_simulated:
        lab.plot(gamma, rates_simulated[channel], 'k+', label="simulated " + str(channel), linewidth=2, markeredgewidth=2)
    lab.xlabel('energy [EeV]')
    lab.ylabel('rate [1/Mpc]')
    lab.legend(loc='lower right')
    lab.grid()
    lab.semilogx()
    lab.savefig('PhotoDisintegration_' + str(id) + '_' + str(channel) +'.png', bbox_inches='tight')
    lab.close()

interaction = mpc.PhotoDisintegration()
table = lab.genfromtxt(mpc.getDataPath('/PhotoDisintegration/pd_table.txt'))
ids = extract_ids(table)

i = 0
count = len(ids)
for id in ids:
    i += 1
    print "Id: ", id, "(", i, "/", count, ")"
    rates_simulated = get_rates_per_channel(interaction, id)
    rates_data = get_data_rates_per_channel(table, id)
    channels = set.union(set(rates_simulated.keys()), set(rates_data.keys()))
    for channel in channels:
        plot_channel(rates_simulated, rates_data, id, channel)
