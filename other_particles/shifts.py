reload(trackage)

big_loop = TL4.loops['u402'].big_loop
other_loop = new_looper

WP = 853807
vz1 = other_loop.tr.p([WP],'z')
vz2 = big_loop.tr.p([WP],'z')
b = trackage.shift_4(nar([vz1]))
c = trackage.shift_4(nar([vz2]))
print(vz1)
print(vz2)
print(vz1-vz2)
