clear all
close all

[Y, FS] = audioread('voce.wav');

%NBITS % quantizazione lineare con elevato numero di bit
lsig = length(Y);
t1  = [0: 1/FS: ((lsig * 1/FS) - 1/FS)]';
%plot(t1, Y);


Nint = 65536;
Graph = 0;
[xi, pdf_x, mean_sig, square_sig, variance_sig] = pdf_estim(Y, Nint, Graph);
%title('pdf segnale campionato e quantizzato a 16 bit');

Z = cumtrapz(pdf_x);

nbits = 8
nbitsfrac = 7
[FIXPTVALUE, DECVALUE, OUTERR, SATFLAG, UNDERFLAG, DIMOP] = FpQuantize(Y, -nbits, nbitsfrac, 'round');

Q = 2^(-nbitsfrac)

Nint = 65536;
Graph = 0;
[xi_2, pdf_n, mean_sig, square_err, variance_sig] = pdf_estim(OUTERR, Nint, Graph);
if Graph == 1 
title('pdf rumore di quantizzazione');
end

%Rapporto SQNR
Psig = square_sig;
Msig = mean_sig;
Perr = square_err;
SQNR_DB_nonComp = 10*log10(Psig / Perr)


NBITS = 16;
NBIT_final = 8;
Nfp_final  = 2^(-(NBIT_final-1));

Nfp_start = 2^(-(NBITS-1));
xstart    = -(2^(NBITS-1)) * Nfp_start;
xstop     = +(2^(NBITS-1) - 1) * Nfp_start;
xvalues   = [xstart: Nfp_start: xstop]'; % possibili valori dell'ingresso xi
xaddress  = [0: 1: (2^NBITS - 1)]';      % indirizzi della tabella look-up

% check look-up table lenght
lut_l = length(xvalues)

% inizialize LUT

%K = +(2^(NBIT_final-1) - 1) * Nfp_final %k in modo tale che il segnale sia quantificabile
Z = Z - Z(length(Z))/2;
Z = Z/(Z(length(Z)));
gz = Z';

%figure
Nint = 100;
Graph = 0;
[xi_3, pdf_n_2, mean_sig_2, square_sig_2, variance_sig_2] = pdf_estim(gz, Nint, Graph);
%pause

% figure
% plot(pdf_x)
% pause
gx = griddedInterpolant(xi, gz);

Nbits = NBIT_final;
nbitsfrac = Nbits - 1; %tutti i bit sulla parte frazionaria
[gx_quant, DECVALUE, OUTERR_LUT, SATFLAG, UNDERFLAG, DIMOP] = FpQuantize(gx(xvalues), -Nbits, nbitsfrac, 'round');

Nfp = Nfp_start;
matlab_offset = 1; %in matlab vectors components start from index 1
Y_addr = round(Y / Nfp) + 2^(NBITS-1); % Y_addr = Equivalente intero + 
                                       % offeset per conversione signed-unsigned
                                       %il round è l'arrotondamento
Y_addr = Y_addr + matlab_offset;

% companding via LUT

Ycomp_LUT = gx_quant(Y_addr);

figure
Nint = 64;
Graph = 1;
[xi_f, pdf_n_f, mean_sig_f, square_sig_f, variance_sig_f] = pdf_estim(Ycomp_LUT, Nint, Graph);
title('pdf segnale campionato e quantizzato a 8 bit - dopo companding via LUT');

Nint = 64;
Graph = 0;
[xi_rf, pdf_n_rf, mean_sig_rf, square_sig_rf, variance_sig_rf] = pdf_estim(OUTERR_LUT, Nint, Graph);

%Rapporto SQNR
Psig = square_sig_f;
Perr = square_sig_rf;
SQNR_DB_Comp = 10*log10(Psig / Perr)
