



APPENDIX W






SOURCE CODE LISTING OF PERSONAL COMPUTER PROGRAM FOR 5 pH POINT TITRATION METHOD



This appendix contains the listing of the source code for the source code file: 5ppt.pas. This file is available on the floppy disk whi·ch is attached to the inside of the thesis cover. The program is coded using Turbo Pascal Ver 4.0. It allows the calculation of H2C03*alkalinity, SCFA (as AT) and systematic pH measurement error from the data collected in the 5 pH point titration procedure. (For 5 pH point titration procedure, see Appendix V).

program atct;

(5 pH POINT TITRATION PROGRAM}

{This program calculates H2C03-alkalinity and short-chain fatty acid concentrations in aqueous solutions containing the carbonate, SCFA, ammoni'-'ll and phosphate subsystems. The total species concentrations of the ammoni'-'ll and phosphate subsystems need to be entered into the c ter program. The sample is titrated from its initial pH to four further pH values measured in the pH region: 6.7, 5.9, 5.2, 4.3. The program provides an estimate of a possible systematic pH measurement error in the event that the carbonate system dominates
over the remaining weak acid/base subsystems.}

{PROGRAM START:}


uses crt;

type phstr= stringt50J;

var pHO, pH1, pH2, pH3, pH4, delph, pHcorr, A1, B1, A2, B2,
HCt1, HCt2, Ct2, Ct1, Ctcomp, HAt1, HAt2, At1, At2, CtAtratio, dil, vx1, vx2, vx3, vx4, vxfi, vxs, ca, temp, ktemp, TDS, logf1, logf2, pK1, pK11, pK2, pK22, pKn, pKnn, Nt,
pKp, pKpp, Pt, pKa, pKaa, At,
HH2C03alksam, H2C03alksam, Hcarb01, Mcarb02, counter, ionstr, nuiber, vsdil, vsundil :real;
again, move,intro, x, y, a, b, c, d: char; labels: array[1..181 of phstr;
value: array[1..18l of real; 1,11,posy: Integer;



function log(x:real):real; begin log:=ln(x)/2.302585093; end;
function tento(y:real):real; begin
 tento:=exp(y1'ln(10)); end;
function nue (TDS,dil:real):real; begin
 Ill.le:= 0.000025* (TOS/dil-20); end;
function logf (mue,ktemp:real):real; begin
logf:= - 1.825 • 1000000 • exp(·1.5*ln(78.3*ktemp)) • (sqrt(111Je)/(1+­ sqrt(111Je)) • 0.3*111.le);
end;
function perH2C03 (ph:real):real; begin

perH2CO3:= 1/(1 + tento(·pk11)/tento(·ph) + tento(·pk11) * tento (·pk22)/tento(-2*ph));
end;
function perHCO3 (ph:real):real;
begin
perHCO3:= 1/(tento(·ph)/tento(·pk11) + 1 + tento(·pk22)/tento (·ph));
end;
function perCO3 (ph:real):real; begin
perC03:=1/(tento(·2*ph)/(tento(·pk11)*tento(·pk22)) + tento C·ph)/tentoC·pk22) + 1);
end;
function per(ph,pkk:real):real; begin
 per:=tento(·pkk)/(tento(·pkk)+tento(·ph)); end;
function dH2CO3alk(pHf,pHs: real): real; begin
 dH2CO3alk:= perHC03(pHf)·perHCO3(pHs) + 2 * (perCO3(pHf)·perC03(pHs)); end;
function dHAcalk(pHf,pHs: real): real; begin
 dHAcalk:= per(pHf,pkaa) • per(pHs,pkaa); end;
function MH20 (vxfi,vxs,pHfi,pHs: real):real; begin
MH2O:= (vsdil+vxs)*tento(·pHs)/tento(logf1) • (vsdil+vxfi)*tento(·pHfi)/ tento(logf1) + Cvsdil+vxfi)*tento(CpHfi-14)) • (vsdil+vxs)* tento((p s-14));
end;
function MNH3(pHf,pHs:real):real; begin
MNH3:= Nt/(14OOO*dil) * vsdil * (per(pHf,pKnn) • per(pHs,pKnn));
end;
function MHP04CpHf,pHs:real):real; begin
   MHP04:= Pt/(31OOO*dil) • vsdil * (per(pHf,pKpp) • per(pHs,pKpp)) end;



procedure mark; begin
textbackground(white); textcolor(black);
end;
procedure urmark;
begin textbackground(black); textcolor(whi te);
end;
procedure box(x1,y1,x2,y2,a:integer); var
 !!:integer; begin
a:= x2-x1·1;
FOR 11:=1 to a DO	{draw the horizontal lines}

begin GOTOXY(x1+1,y1);
writeln(chr(196)); GOTOXY(x1+1,y2);
writeln(chr(196));
x1:=x1+1;
end;
x1:=x1·a; a:=y2-y1·1;
FOR 11:=1 to a DO
begin GOTOXY(x1,y1+1);
writeln(chr(179)); GOTOXY(x2,y1+1);
writeln(chrC179)); y1:=y1+1;
  end; y1:=y1·a;
GOTOXY(x1,y1)j
writeln(chrC218)); GOTOXY(x1,y2);
writeln(chr(192)); GOTOXY(x2,y1)j
writeln(chr(191)); GOTOXY(x2,y2);
   writeln(chr(217)); end; {procedure box}
procedure introscreen; begin
clrscr; mark;








{draw the vertical lines}






{draw the corners}

box(10,3,70,20,1); urmark; box(24,6,56,8,1);box(22,5,58,9,1); box(26,10,55,16,1);	-
GOTOXY(25,7); highvideo;
write(1FIVE pH POINT TITRATION METHOO');normvideo; GOTOXY(28,11);
write(1FOR DETERMINATION OF:1); GOTOXY(28,13);
write('(1) SHORT-CHAIN FATTY ACIDS'); GOTOXY(28,15);
write('C2) H2C03*ALKALINITY1); GOTOXY(11,16);
GOTOXY(12,18);
write('Copyright: UCT'); box(10,22,41,24,1); GOTOXY(12,23);
write('Press any letter to continue'); intro:= readkey;
end; {procedure introscreen}
procedure initlabel; begin
labelsC1J:='pHo (initial pH)•••••••••••••••••; labelsC2J:='pH1 (after adding Vx1)	•;
labels[3J:='pH2 (after adding Vx2)	•;
labels[4]:='pH3 (after adding Vx3)	•;
labelsC5J:=1pH4 (after adding Vx4)	•;
labels[6J:=1	Vx1 (ml)	•;
labelsCn:=1	Vx2 (ml)	•;
labelsC8J:='	Vx3 (ml)	•;
labels[9J:=•	Vx4 (ml)	•;
labelsC10J:='Normality of titrant (mol/l)	•;
labelsC11J:=1Sample size: undiluted (ml)	•;
labelsC12J :='Sample size: diluted (ml)	•;
labelsC13J :='Terrperature (Celsius)	•;
labelsC14J:='TDS (mg/l)	•;
labels[15J:=1Specific Conductivity (mS/m)	•;
labels[16J:='lnorganic Nitrogen (mgN/l)	•;

"· 5
labelsc1n:=1lnorganic Phosphorus (mgP/l)	';
end; {procedure initlabel}
procedure default_values; begin
value[1]:= 7.36;
value[2]:= 6.75;
value[31:= 5.95;
value[4]:= 5.18;
value[5J:= 4.29;
value[6]:= 1.06;
value[7]:= 3.50;
value[8]:= 4.84;
value[91:= 5.40;
value[10]:= 0.0728;
value[111 := 10;
value[12]:= 50;
value[13]:= 21;
value[14]:= 3300;
value[15J:= 488; value[16]:= O; valueC1n:= O;
end; {procedure default_values}
procedure allocate_data;
begin
pHO:= value[1l; pH1:= value[2]; pH2:= value[3]; pH3:= value[4]; pH4:= value[5J; vx1:= value[6]; vx2:= valueen; vx3:= value[8]; vx4:= value[9J; ca:= value[10J; vsundil:= value[11J; vsdil:= value[12]; tefl1):= value[13]; TOS:= value[14J; ionstr:= value[15]; Nt:= value[16J;
Pt:= valueC1n;
end; {allocate_data}
procedure write_value; begin GOTOXYC37,2+posy);
case posy of
1,2,3,4,5,6,7,8,9: begin
writeln(value[posy]:1:2);
end;
10: begin
 writeln(value[posy]:1:4); end;
11,12,13,14,15,16,17,18: begin
writeln(value[posy]:1:0);
end; else end;
end; {procedure write_value}
procedure screen; begin
mark; box(50,2,79,8,1); urmark;
box(50,10,79,19,1);
boxC50,10,59,19,1);
box(50,10,79,13,1); GOTOXY(52,4);

11.6

write(' 5 pH POINT TITRATION '); GOTOXYC52,6);
write(' TITRATION INPUT DATA'); GOTOXY(52,12);
writeln(' KEY	FUNCTION'); GOTOXY(52,14);
writeln('	1,CHR(24),CHR(25),' SELECT PARAMETER');
GOTOXYC52,15);
writeln('<enter> ERASE VALUE'); GOTOXY(52,16);
writeln('<enter> INSERT NEIi VALUE'); GOTOXY(52,17);
writeln(' C	CALCULATION'); GOTOXY(52,18);
 writeln(' Q	QUIT'); end; {procedure sreeen)
procedure restore;
begin clrscr; screen;
mark; box(1,1,47,22,1); unmark;
FOR I:= 1 to 17 Do
begin GOTOXY(3,2+I);
writeln(labels[ll ); GOTOXY(37,2+1);
case I of
1,2,3,4,5,6,7,8,9: begin
writeln(va Lue[ll:1:2);
end;
1D: begin
writeln(value[l]:1:4);
end;
11,12,13,14,15,16,17,18: begin
writeln(value [I]:1:0);
end;
else end;
end;  {FOR)
GOTOXY(3,2+posy);
mark; write(labels[posy]); write value;
 unmark; end;
procedure display;
begin clrscr; screen;
mark; box(1,1,47,21,1); unnark;
FOR I:= 1 to 17 Do
begin
if 1=1 then mark;
GOTOXY(3,2+1);
writeln(labels[ll); GOTOXY(37,2+1);
case I of
1,2,3,4,5,6,7,8,9: begin
writeln(value[I]:1:2);
end;
10: begin
writeln(value(Il:1:4);
end;
11,12,13,14,15,16,17,18: begin
writeln(value[Il:1:0);
end;
else

end;
if 1=1 then unmark;
 end; posy:=1; REPEAT
repeat move:= readkey until (upcase(move) in [1C1,1Q1]) or (ORD(move)
in [72,80,13]);
case 0RD(move) of 80: begin
G0T0XY(3,2+posy);
writeln(labels[posyl); write value; posy:;posy+1;
if posy=18 then posy:=1; G0T0XY(3,2+posy);
 mark; writeln(labels[posy]); write_value; unmark; end;
72: begin
G0T0XY(3,2+posy);
writeln(labels[posyJ); write value; posy:;posy-1;
if posy=O then posy:=17; G0T0XY(3,2+posy);
 mark; writeln(labels[posyJ); write_value; unmark; end;
13: begin
repeat G0T0XY(37,2+posy);
mark;
writeln('	'); GOT0XY(37,2+posy);
 ($!·} read(nuit>er) ($1+} until ioresult = O; value[posy]:= nunber; G0T0XY(37,2+posy);
unmark;
writeln('	'); write value;
case posy of 14: begin
value[posy+1J:= O.1488*(value[posyJ-20);
if value[posy+1J < 0 then value[posy+1J:= 0; GOTOXY(37,posy+3);
write('	'); G0T0XY(37,posy+3);
 write(value[posy+1]:1:0); end;
15: begin
value[posy-1]:= 6.n•value[posyJ+20; G0T0XY(37,posy+1);
write('	'); GOT0XY(37,posy+1);
 write(value[posy-1]:1:0); end;



else
end;
 
else end; end;

until upcase(move) in c•c1,1Q1];
clrscr; writeln('Calculating, please wait');
end; (procedure display}

procedure p(; begin
kt := 273 + t ;

if TDS < 20 then TDS:= 21; logf1:= logf(111Je(TDS,dil),kt ); logf2:= 4*logf1;
pK1:= ·1*(-356.3094 - 0.06091964*kt + 21834.37/kt + 126.8339
     *log(kt ) - 1684915/(kt *kt )); pK11:= pK1 + logf1;
pK2:= ·1*(·107.8871 - 0.03252849*kt + 5151.79/kt + 38.92561
      *log(kt ) - 563713.9/(kt *kt >>; pK22:= pK2 - logf1 + logf2;
pKa:= 1170.5/kt - 3.165 + 0.0134 * kt ;
pKaa:= pKa + logf1;	.
pKn:= 2835.8/kt -0.6322 + 0.00123 * kt ;
pKnn:= pKn + logf1 ;
pKp:= 1979.5/kt - 5.3541 + 0.01984 * kt ;
 pKpp:= pKp - logf1 + logf2 end; {procedure pK}
procedure deltapH; begin

 pHO:= pHO + pHcorr; pH3:= pH3 + pHcorr; end;

pH1:= pH1 + pHcorr; pH4:= pH4 + pHcorr;

pH2:= pH2 + pHcorr;

procedure atctcalculation; begin
A1:= (vx2-vx1)*ca • MH20(vx1,vx2,pH1,pH2) - MNH3(pH1,pH2) - MHP04(pH1,pH2) + dHAcalk(pH1,pH2)/dHAcalk(pH3,pH4) * (MH20(vx3,vx4,pH3,pH4) + MNH3(pH3,pH4) + MHP04(pH3,pH4) -
(vx4-vx3)*ca);
B1:= dH2C03alk(pH1,pH2) • dHAcalk(pH1,pH2)/dHAcalk(pH3,pH4) * dH2C03alk(pH3,pH4);
A2:= (vx4-vx1)*ca - MH20(vx1,vx4,pH1,pH4) - MNH3(pH1,pH4) -
 MHP04(pH1,pH4) + dHAcalk(pH1,pH4)/dHAcalk(pH3,pH4) * (MH20(vx3,vx4,pH3,pH4) + MNH3(pH3,pH4) + MHP04(pH3,pH4) •
(vx4-vx3)*ca);
B2:= dH2C03alk(pH1,pH4) - dHAcalk(pH1,pH4)/dHAcalk(pH3,pH4) * dH2C03alk(pH3,pH4);

Ct1:= (A1/B1)/vsdil * 50000 * dil; Ct2:= (A2/B2)/vsdil * 50000 * dil; MCt1:= A1/B1;
Ctc°"":= Ct1 • Ct2;
end; {procedure atctcalculation}
procedure atct1; begin
dil:= vsdil/vsundil;

{us,ing pH: 1-2; 3-4}
{using pH: 1-4; 3-4}

pK; delpH:= O; pHcorr:=0; atctcalculation; counter:= O;
MAt1:= 1/dHAcalk(pH3,pH4) * ((vx4·vx3)*ca·MCt1*dH2C03alk(pH3,pH4) • MNH3(pH3,pH4) - MHP04(pH3,pH4) • MH20(vx3,vx4,pH3,pH4));
At1:= MAt1/vsdil * 60000 * dil; CtAtratio:= At1/Ct1;
if etc°""= 0 then x:= 1d1; if Ctc°"" > 0 then x:= •a•; if etc°""< 0 then x:= 1b1;
if CtAtratio > 0.5 then x:= •c•;
case x of
•a•: begin
pHcorr:= -0.01; repeat
delpH:= delpH • 0.01;
deltapH; counter:= counter+ 1; atctcalculation;
until (Ctc°"" < 0) or (counter> 19);
end;
1b1:  begin
pHcorr:= 0.01; repeat





else
end;

delpH:= delpH + 0.01;
deltapH; counter:= counter+ 1; atctcalculation;
until (Ctcomp > 0) or (counter> 19);
end;

MAt1:= 1/dHAcalk(pH3,pH4) * ((vx4-vx3)*ca-MCt1*dH2C03alk(pH3,pH4) - MNH3(pH3,pH4) • MHP04(pH3,pH4) - MH20(vx3,vx4,pH3,pH4));
At1:= MAt1/vsdil * 60000 * dil;
MAt2:= 1/dHAcalk(pH1,pH4) * ((vx4-vx1)*ca·MCt1*dH2C03alk(pH1,pH4) - MNH3(pH1,pH4) • MHP04(pH1,pH4) • MH20(vx1,vx4,pH1,pH4));
At2:= MAt2/vsdil * 60000 * dil;
MH2C03alksam:= MCt1 * (perHC03(pH0) + 2 * perC03CpH0));
H2C03alksam:= MH2C03alksam / vsdil * dil * 50000 + (tento(CpH0-14)) - tento(-pH0)/tento(logf1)) * 50000;
end; {procedure atct1}
procedure output; begin
clrscr;
mark; box(10,2,60,4,1);unmark; box(2,7,73,10,1); box(49,7,73,10,1);
box(2,11,73,14,1); box(49,11,73,14,1);
box(2,15,73,18,1); box(49,15,73,18,1); G0T0XY(26,3);
highvideo; write('OUTPUT DATA'); normvideo; G0T0XY(5,8);
write('H2C03*alkalinity (undiluted saq,le)'); G0TOXY(58,8);
highvideo; write(H2C03alksam:3:0); normvideo;
G0T0XY(5,9);	1
write('(mg/l as CaC03) ); G0T0XY(5,12);
write('Short-chain fatty acids (undiluted saq,le)'); G0TOXY(58,12);
if At1 > 0 then begin
 highvideo; write(At1:3:0); normvideo; end
else begin
 highvideo; write('0 ');normvideo; end;
G0TOXY(5,13);
write('(mg/l as acetic acid)'); G0TOXY(5,17);
write(1Systematic pH error•); G0TOXY(58,17);
highvideo; write(delpH:1:2); normvideo; if counter= 20 then
begin
G0TOXYC5,16);
write('The titration data indicate a systematic'); G0TOXYC5,17);
writeC'pH error> 0.2; Check pH probe calibration'); G0TOXYC58,17);
write('	1 );
end;
if x = •c• then begin G0TOXYC5,16);
write('Correction for systematic pH error'); G0TOXY(5,17);
write('is not possible for this titration'); G0T0XY(60,17); •

 write('	'); end;
end; {procedure output}
procedure exitbox;
begin
repeat;
box(S,20,60,24,1); box(49,20,60,22,1);
box(S,22,60,24,1); box(49,22,60,24,1); GOTOXY(7,21 );
write('Do you wish to do a further calculation 71); GOTOXY(53,21);
highvideo; write('Y'); normvideo; write('/'); highvideo; write ('N');
normvideo;
GOTOXY(7,23);
write('Do you wish to quit the program?'); GOTOXY(54,23);
 highvideo; write('Q'); normvideo; again:=readkey; again:=upcase(again); until again in [1Y1,1N1,1Q1J;
clrscr;
end; {procedure exitbox}


begin	C of main program} introscreen;
again:= 'Y'; default_values;
repeat; clrscr; initlabel; display; if upcase(move) in ['C'l then begin
 allocate_data; atct1; output; exitbox; end
else begin
clrscr; again:= 'N';
-end;
until agein in c•N 1 , 1 Q•J ;
end. {of main program}

{PROGRAM END.}


















