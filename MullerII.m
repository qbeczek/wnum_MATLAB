%*ustawienia poczatkowe*
close all
clear all

deltax = 10^-3;
w_const = [-72090, 1458, -4536, -1152, -14, 14, 1];
MAX_ITER = 100000;
w_cur = w_const;

%*muller wersja II, pierwiastki szukane: 1+3j, 1-3j*

x_i= 45 + 120*i;
x_dec1= 30 + 90*i;
x_dec2= 20 + 60*i;
iter = 0; %pocz�tek p�tli iteracyjnej metody mullera wersji II
while (iter<MAX_ITER)
    iter = iter + 1;
    UP = [-(x_i - x_dec1)^2 (x_i - x_dec1);-(x_i - x_dec2)^2 (x_i - x_dec2)];
    DOWN = [w(w_cur, x_i)-w(w_cur, x_dec1); w(w_cur, x_i)-w(w_cur, x_dec2)];
    mat = UP\DOWN;              %lewe dzielenie macierzy
    ai = mat(1);
    bi = mat(2);
    ci = w(w_cur, x_i);
    x_inc1 = x_i - (2*ci)/(bi+sign(bi)*sqrt(bi^2 - 4*ai*ci));
    x_dec2 = x_dec1;
    x_dec1 = x_i;
 
    if (abs(x_inc1-x_i)<deltax) break;
    end
    
    x_i = x_inc1;
end
iter
x1 = x_inc1
x2 = x_inc1'
w_cur = def_kw(w_cur, x_i);


%*Muller wersja I, pierwiastki do znalezienia: -8 + 5j, -8 - 5j*

clear x_i x_inc1 ai bi ci iter
x_i = -34 + 25*i;
wbx_diff1 = pochodna(w_cur, 1);
wbx_diff2 = pochodna(w_cur, 2);
iter = 0; %pocz�tek p�tli iteracyjnej metody mullera wersji I
while (iter<MAX_ITER)
    iter = iter + 1;
    ai =@(x_i) ((1./2).*w(wbx_diff2, x_i));
    bi =@(x_i) (w(wbx_diff1, x_i));
    ci =@(x_i) (w(w_cur, x_i));
    x_inc1 = x_i - (2*ci(x_i))/(bi(x_i)+sign(bi(x_i))*sqrt(bi(x_i)^2 - 4*ai(x_i)*ci(x_i)));
    
    if (abs(x_inc1-x_i)<deltax) break;
    end
   
    x_i = x_inc1;
end
iter
x3 = x_inc1
x4 = x_inc1'
w_cur = def_kw(w_cur, x_i);

%*Metoda siecznych, pierwiastek do znalezienia: 9
clear iter x_i x_inc1
x_i = 9;
x_dec1 = 8;
iter = 0; %pocz�tek p�tli iteracyjnej metody siecznych
while (iter<MAX_ITER)
    iter = iter + 1;
    x_inc1 = x_i - (x_i-x_dec1)/(w(w_cur, x_i)-w(w_cur, x_dec1))*w(w_cur, x_i);
    x_dec1 = x_i;
    
    if (abs(x_inc1-x_i)<deltax) break;
    end
    x_i = x_inc1;
end
iter
x5 = x_inc1
w_cur = def_lin(w_cur, x_i)

%*metoda stycznych, pierwiastek do znalezienia: -9
clear x_i iter x_inc1
wbx_diff = pochodna(w_cur, 1);
x_i = -30;
iter = 0; %%pocz�tek p�tli iteracyjnej metody stycznych Newtona
while (iter<MAX_ITER)
    iter = iter + 1;
    x_inc1 = x_i - w(w_cur, x_i)/w(wbx_diff, x_i);
    if (abs(x_inc1-x_i)<deltax) break;
    end
    x_i = x_inc1;
end
iter
x6 = x_inc1

%FUNKCJE

%*funckja obliczaj�ca pochodn� wielomianu dowolnego stopnia*
function w_diff = pochodna(w_old, i)
w_tmp = w_old;
for k = 1 : i
    w_tmp = polyder(w_tmp); %Funkcja polyder(p) oblicza wektor wsp�czynnik�w wielomianu b�d�cego
                            %pochodn� wielomianu reprezentowanego przez p.
end
w_diff = w_tmp;
end

%*deflacja kwadratowa*
function w_new = def_kw(w_old, x)
p = 2*real(x);
r = -(abs(x)*abs(x));
n = length(w_old);     %stopien wielomianu +1
w_new(n-2) = w_old(n); %deflacja obniza stopien o 2
w_new(n-3) = w_old(n-1)+p*w_new(n-2);
   for n=(n-4):-1:1
    w_new(n) = w_old(n+2) + p*w_new(n+1)+r*w_new(n+2);
   end
end

%*deflacja liniowa*
function w_new = def_lin(w_old, x)
w_new(length(w_old) - 1) = w_old(length(w_old))
    for n = (length(w_new)-1):-1:1
        w_new(n) = w_old(n+1) + x*w_old(n+1);
    end
end


%*funkcja wielomianowa*
function w_new = w(w_old, x)
k = length(w_old);
w_new = 0;
    for n=0 : (k-1)
        w_new = w_new + w_old(n+1)*x^n;
    end
end
