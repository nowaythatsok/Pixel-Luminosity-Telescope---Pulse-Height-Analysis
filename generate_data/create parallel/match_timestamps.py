

fill_ID_to_timestamp_dict = {}
f_fill_ID_to_timestamp = open("fill_ID_to_timestamp.dat", "r")
for line in f_fill_ID_to_timestamp:

    split_line = line.split("\t")
    fill_ID_to_timestamp_dict[split_line[0]] =  split_line[1].split(".")

f_fill_ID_to_timestamp.close()



timestamps_available = []
f_timestamps = open("timestamps_available.dat", "r")
for line in f_timestamps:

    split_line = line.split("_")[1].split(".")
    timestamps_available += [[split_line[0], split_line[1]]]

f_timestamps.close()

problematic = []
fill_IDs = sorted(fill_ID_to_timestamp_dict.keys())
for fill_ID in fill_IDs:
    ts1     = fill_ID_to_timestamp_dict[fill_ID]
    date1   = ts1[0]
    time1   = ts1[1]

    D = []
    lTs= []

    for ts2 in timestamps_available:
        if date1 == ts2[0]:
            D+=[abs(int(time1)-int(ts2[1]))]
            lTs+=[ts2]

    if len(D)>0:
        print 'VVV_peters_pulse_height_controller ' '/localdata/2017/SLINK/Slink_'+".".join(str(_) for _ in lTs[D.index(min(D))])+'.dat   ./GainCal/GainCalFits_20170518.143905.dat  ' + fill_ID
    else:    
        #print "warning", fill_ID
        problematic += [fill_ID]

#parallel.sh




