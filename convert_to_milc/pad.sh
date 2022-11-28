for i in {0..9}
   do dd if=/dev/zero bs=1 count=4096 >> l2464f211b600m0102m0509m635a.000100$i
done

for i in {10..99}
   do dd if=/dev/zero bs=1 count=4096 >> l2464f211b600m0102m0509m635a.00010$i
done

for i in {100..143}
   do dd if=/dev/zero bs=1 count=4096 >> l2464f211b600m0102m0509m635a.0001$i
done
