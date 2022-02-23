function y = chi2(xi)
%% 
%%  Calculation of VCKMs and chi^2 
%%
Au = xi(:,1);
Ad = xi(:,2);
EEu = xi(:,3);
EEd = xi(:,4);
phi1 = xi(:,5);
phi2 = xi(:,6);

        y=chi2CKM23(Au, Ad, EEu, EEd, phi1, phi2)+...
            chi2CKM12(Au, Ad, EEu, EEd, phi1, phi2)+...
            chi2CKM13(Au, Ad, EEu, EEd, phi1, phi2)+...
            chi2JJ(Au, Ad, EEu, EEd, phi1, phi2);
end
%%/*Se definen los parámetros de la siguiente manera: */
function y = BB(A,EE,lambda1,lambda2,lambda3)
        y=sqrt(((A-lambda1)*(A-lambda2)*(lambda3-A))/(A-EE));
end

function y = DD(A,EE,lambda1,lambda2,lambda3)
        y=sqrt(((EE-lambda1)*(lambda2-EE)*(lambda3-EE))/(A-EE));
end

%%/* los elementos diagonales */
function y = O11(A,EE,lambda1,lambda2,lambda3)
        y=1/(sqrt(1+((EE-lambda1)*(A-EE)/((lambda2-EE)*(lambda3-EE)))+((EE-lambda1)*(A-lambda2)*(lambda3-A)/((lambda2-EE)*(lambda3-EE)*(A-lambda1)))));
end
function y = O22(A,EE,lambda1,lambda2,lambda3)
        y=1/(sqrt(1+((EE-lambda1)*(lambda3-EE)/((lambda2-EE)*(A-EE)))+((A-lambda1)*(lambda3-A)/((A-EE)*(A-lambda2)))));
end

function y = O33(A,EE,lambda1,lambda2,lambda3)
        y=1/(sqrt(1+((lambda3-A)*(A-EE)/((A-lambda2)*(A-lambda1)))+((EE-lambda1)*(lambda2-EE)*(lambda3-A)/((A-lambda2)*(lambda3-EE)*(A-lambda1)))));
end

%%/*los elementos fuera de la diagonal, donde todos dependen de A, EE, lambda1, lambda2, lambda3*/
function y = O12(A,EE,lambda1,lambda2,lambda3)
        y=abs(DD(A, EE, lambda1, lambda2, lambda3))/(lambda2-EE)*O22(A, EE, lambda1, lambda2, lambda3);
end

function y = O21(A,EE,lambda1,lambda2,lambda3)
        y=(-(EE-lambda1)/(abs(DD(A, EE, lambda1, lambda2, lambda3))))*O11(A, EE, lambda1, lambda2, lambda3);
end

function y = O32(A,EE,lambda1,lambda2,lambda3)
        y=(-(abs(BB(A, EE, lambda1, lambda2, lambda3)/(A-lambda2))))*O22(A, EE, lambda1, lambda2, lambda3);
end

function y = O23(A,EE,lambda1,lambda2,lambda3)
        y=((lambda3-A)/abs(BB(A, EE, lambda1, lambda2, lambda3)))*O33(A, EE, lambda1, lambda2, lambda3);
end

function y = O13(A,EE,lambda1,lambda2,lambda3)
        y=abs(DD(A, EE, lambda1, lambda2, lambda3)/(lambda3-EE))*((lambda3-A)/(abs(BB(A, EE, lambda1, lambda2, lambda3))))*O33(A, EE, lambda1, lambda2, lambda3);
end

function y = O31(A,EE,lambda1,lambda2,lambda3)
        y=(abs(BB(A, EE, lambda1, lambda2, lambda3))/(A-lambda1))*((EE-lambda1)/abs(DD(A, EE, lambda1, lambda2, lambda3)))*O11(A, EE, lambda1, lambda2, lambda3);
end

%%/* Los elementos de la VCKM */
function y = CKM11(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O11(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb)) +...
            (O21(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))) +...
            (O31(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

function y = CKM22(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O12(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+...
            O22(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            O32(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2));
end

function y = CKM33(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O13(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+...
            O23(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            O33(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2));
end

function y = CKM12(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O11(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+...
            O21(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            O31(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2));
end

function y = CKM21(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O12(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb))+...
            O22(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            (O32(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

function y = CKM13(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O11(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+...
            O21(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            (O31(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

function y = CKM31(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O13(Au, EEu, mu, mc, mt)*O11(Ad, EEd, md, ms, mb))+...
            O23(Au, EEu, mu, mc, mt)*O21(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            (O33(Au, EEu, mu, mc, mt)*O31(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

function y = CKM32(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O13(Au, EEu, mu, mc, mt)*O12(Ad, EEd, md, ms, mb))+...
            O23(Au, EEu, mu, mc, mt)*O22(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            (O33(Au, EEu, mu, mc, mt)*O32(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

function y = CKM23(Au,Ad,EEu,EEd,phi1,phi2)
global mu;
global md;
global ms;
global mc;
global mb;
global mt;
        y=(O12(Au, EEu, mu, mc, mt)*O13(Ad, EEd, md, ms, mb))+...
            O22(Au, EEu, mu, mc, mt)*O23(Ad, EEd, md, ms, mb)*exp(complex(0,phi1))+...
            (O32(Au, EEu, mu, mc, mt)*O33(Ad, EEd, md, ms, mb)*exp(complex(0,phi1+phi2)));
end

%%/*El invariante de Jarlskog*/
function y = JJ(Au,Ad,EEu,EEd,phi1,phi2)
        y=imag(CKM12(Au, Ad, EEu, EEd, phi1, phi2)*CKM23(Au, Ad, EEu, EEd, phi1, phi2)*...
            conj(CKM13(Au, Ad, EEu, EEd, phi1, phi2))*conj(CKM22(Au, Ad, EEu, EEd, phi1, phi2)));
end

function y = chi2CKM23(Au,Ad,EEu,EEd,phi1,phi2)

global vcb_th;
global sigma_vcb;

        y=(vcb_th -abs(CKM23(Au, Ad, EEu, EEd, phi1, phi2)))^2/(sigma_vcb)^2;
end

function y = chi2CKM12(Au,Ad,EEu,EEd,phi1,phi2)

global vus_th;
global sigma_vus;

        y=(vus_th -abs(CKM12(Au, Ad, EEu, EEd, phi1, phi2)))^2/(sigma_vus)^2;
end

function y = chi2CKM13(Au,Ad,EEu,EEd,phi1,phi2)

global vub_th; 
global sigma_vub;

        y=(vub_th -abs(CKM13(Au, Ad, EEu, EEd, phi1, phi2)))^2/(sigma_vub)^2;
end

function y = chi2JJ(Au,Ad,EEu,EEd,phi1,phi2)

global jarslkog_th;
global sigma_jars;

        y=(jarslkog_th -JJ(Au, Ad, EEu, EEd, phi1, phi2))^2/(sigma_jars)^2;
end
