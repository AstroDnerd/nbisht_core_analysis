

core_list=[12,13,14,15,16,17,18,19,21]
Nc = len(core_list)
cores_used = nar(ht1.cores_used)
fff = " [%2d]"*Nc
head = "    " + fff%tuple(core_list)
print(head)
for core_id in core_list:
    overlap = nar(ht1.overlaps[core_id])
    these_cores = cores_used[ overlap>0]
    these_overlaps = overlap[ overlap > 0]
    drr = dict(zip(these_cores,these_overlaps))
    output= "[%2d]"%core_id
    for ccc in core_list:
        if ccc in drr:
            output += " %0.2f"%drr[ccc]
        else:
            output += " ----"
    print(output)
