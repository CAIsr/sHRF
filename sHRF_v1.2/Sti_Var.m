function Sti_new = Sti_Var(STI, Vars)

Sti_new = zeros(size(STI));
tau = round(Vars(1)*10);
Sti_new(:,1) = STI(:,1);

ii =0;
for i = 1:1:(size(STI,1)-tau)
    if STI(i,2) > 0
        jj = floor(ii/3)+1;
        if mod(jj,16) == 0
            kk = 17;
        else
            kk = mod(jj,16)+1;
        end
        Sti_new((i+tau),2)=STI(i,2) * (1+Vars(kk));
        ii = ii + 1;
    end
end
end

