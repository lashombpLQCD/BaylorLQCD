function [thq] = ModLejaComplex(th)


pp = length(th);

    thr(1) = th(1);

    index(1) = 1;

    j = 2;

    for i=2:pp

        if( abs(th(i)) > abs(thr(1)) )

          thr(1) = th(i);

        end

    end

    iatendofloop = i

    while( j <= pp )

        for i=1:pp

            pr(i) = 0;

            for ii=1:j-1

                pr(i) = pr(i) + log(abs(th(i)-thr(ii)));

            end

        end

        thr(j) = th(1);

        prj = pr(1);

        itemp = 1;

        for i=2:pp

            if( pr(i) > prj )

                thr(j) = th(i);

                prj = pr(i);

            end

        end

        j = j + 1;

    end

  

    thq = transpose(thr);