import mpc
import pylab as lab
import sys
import math
# make sure root is startet in Batch mode
sys.argv.append('-b')
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.StartGuiThread = False 
ROOT.gStyle.SetOptFit(1)

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
        if i % 100000 == 0:
            print ".",
        c.setNextStep(10000000 * mpc.Mpc)
        c.clearInteractionStates()
    
        interaction.process(c)
        state = mpc.InteractionState()
        c.getInteractionState("mpc::PhotoDisintegration", state)
        if not h.has_key(state.channel):
            h[state.channel] = []
        v = state.distance / mpc.Mpc
        h[state.channel].append(v)
        
    for key in h.keys():
        l = h[key]
        minV = min(l)
        maxV = max(l)
        hist = ROOT.TH1F(str(id) + str(key) + str(gamma), '', 20, minV, maxV)
        for v in l:
            hist.Fill(v)
        h[key] = hist
        
    return h

def get_rates_per_gamma(interaction, id, gamma):
    r = {}
    h = get_rate_h(interaction, id, gamma, 100000)
    print "create hists"
    for key in h.keys():
        hist = h[key] 
        if hist.GetEffectiveEntries() == 0:
            r[key] = 0
            continue
        f = ROOT.TF1('f1', 'expo')
        hist.Fit(f, "q")
        c = ROOT.TCanvas()
        #file = ROOT.TFile(str(id) + "_" + str(gamma) + "_" + str(key) + ".root", "CREATE OVERWRITE")
        hist.Draw()
        #hist.Write()     
        c.Print('PhotoDisintegration_' + str(id) + '_' + str(key) + '_' + str(lab.log10(gamma)) + '.png')
        c.Close()
        #file.Close()
        slope = -f.GetParameter(1)
        if not (slope > 0 and slope < 100):
            slope = 0.00000001
        r[key] = slope
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
        if iE % 20 != 0:
            continue
        rates_per_gamma = get_rates_per_gamma(interaction, id, gamma[iE])
        for channel in rates_per_gamma.keys():
            if not rates.has_key(channel):
                rates[channel] = lab.zeros(200)
            rates[channel][iE] = rates_per_gamma[channel]
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
        y = list(rates_simulated[channel])
        print y 
        lab.plot(gamma, y, 'k+', label="simulated " + str(channel), linewidth=2, markeredgewidth=2)
    lab.xlabel('energy [EeV]')
    lab.ylabel('rate [1/Mpc]')
    lab.semilogy()
    lab.ylim(1e-9, 100)
    lab.legend(loc='lower right')
    lab.grid()
    lab.semilogx()
    lab.savefig('PhotoDisintegration_' + str(id) + '_' + str(channel) + '.png', bbox_inches='tight')
    lab.close()

interaction = mpc.PhotoDisintegration()
table = lab.genfromtxt(mpc.getDataPath('/PhotoDisintegration/pd_table.txt'))
ids = [1044019000] #extract_ids(table)

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
    break
