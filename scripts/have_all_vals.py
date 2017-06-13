runfiles = "PartDet/"

compare = "Analyses/example/"

files = ["Cuts.in", "DiParticle_info.in", "ElectronTau_info.in", "Electron_info.in", "Gen_info.in", "Jet_info.in", "FatJet_info.in", "MuonElectron_info.in", "MuonTau_info.in", "Muon_info.in", "Run_info.in", "Tau_info.in", "VBFCuts_info.in"]

for onefile in files:
    f = open(compare+onefile, 'r')
    temp = f.read().splitlines()
    mapper = {}
    for line in temp:
        if len(line) == 0:
            continue
        if line[0] == '#' or (line[0] == '/' and line[1] == '/'):
            continue
        splitline = line.split()
        if len(splitline) <= 1:
            continue
        mapper[splitline[0]] = False

    foundnotinrun = False
    f2 = open(runfiles+onefile, 'r')
    temp = f2.read().splitlines()
    for line in temp:
        if len(line) == 0:
            continue
        if line[0] == '#' or (line[0] == '/' and line[1] == '/'):
            continue
        splitline = line.split()
        if len(splitline) <= 1:
            continue
        if splitline[0] not in mapper:
            if not foundnotinrun:
                print onefile
                print "--------------------------------"
                print "In Run files, but not example"
                foundnotinrun = True
            print splitline[0]
        else:
            mapper[splitline[0]] = True
    foundnotinex = False
    for val, found in mapper.items():
        if not found:
            if not foundnotinex and not foundnotinrun:
                print onefile
                print "--------------------------------"
                print "In Example files, but not run"
                foundnotinex = True
            if not foundnotinex:
                print
                print "In Example files, but not run"
            print val

    if foundnotinex or foundnotinrun:
        print

