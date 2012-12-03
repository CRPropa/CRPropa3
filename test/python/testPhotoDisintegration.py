from mpc import *
from pylab import *


def get_rates(candidate, N=5000):
    # perform N interactions
    D = {}
    for i in range(N):
        interaction = InteractionState()
        pdModule.setNextInteraction(candidate, interaction)
        if interaction.channel == 0:
            continue
        D.setdefault(interaction.channel, []).append(interaction.distance)

    # calculte exclusive mean interaction lengths
    for k in D.keys():
        l = array(D[k])
        p = len(l) / float(N) # probability for channel
        l = mean(l) / p / Mpc # exclusive decay length in Mpc
        D[k] = 1 / l # exclusive decay rate in 1/Mpc
    return D

def get_data_rates(pid):
    D = {}
    for i in range(len(table)):
        z, n = int(table[i][0]), int(table[i][1])
        if nucleusId(z+n, z) != pid:
            continue
        channel = table[i][2]
        rates = table[i][3:]
        D[channel] = rates
    return D

def parse_pid(pid):
  return 'Z=%i, A=%i'%((pid//10000)%1000, (pid%10000)/10)

def parse_channel(channel):
    s = '%06d' % channel
    d = list(map(int, s))
    s = 'n, ' * d[0] + 'p, ' * d[1] + 'H$^2$, ' * d[2] + 'H$^3$, ' * d[3] + 'He$^3$, ' * d[4] + 'He$^4$, ' * d[5]
    return s[0:-2]



### nucleus to test
pid = nucleusId(4,2)

candidate = Candidate()
candidate.current.setId(pid)

table = genfromtxt(getDataPath('photodis_CMB.txt'))
pdModule = PhotoDisintegration(CMB)
gamma = logspace(6, 14, 200)
gamma2 = gamma[1::5]

Dsim = {}
for i,g in enumerate(gamma2):
    candidate.current.setLorentzFactor(g)
    D = get_rates(candidate)
    for k in D.keys():
        Dsim.setdefault(k, zeros(len(gamma2)))[i] = D[k]

Ddata = get_data_rates(pid)

# plot
for k in Dsim.keys():
    figure()
    plot(gamma2, Dsim[k], 'k+')
    plot(gamma, Ddata[k], 'b')
    s = 'Nucleus ' + parse_pid(pid) + '\nDisintegration Channel ' + parse_channel(k)
    text(0.1, 0.85, s, transform=gca().transAxes)
    xlabel('Lorentzfactor $\gamma$')
    ylabel('Rate [1/Mpc]')
    ylim(1e-6, 1e3)
    loglog()
    grid()
    savefig('PhotoDisintegration_' + str(pid) + '_' + str(k) + '.png', bbox_inches='tight')
    show()




