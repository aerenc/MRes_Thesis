function X_xi = get_Dummies (xxx,mmm,ttt)

xxx_1  = xxx==1;
xxx_1  = double(xxx_1);
xxx_2  = xxx==2;
xxx_2  = double(xxx_2);
xxx_3  = xxx==3;
xxx_3  = double(xxx_3);
xxx_4  = xxx==4;
xxx_4  = double(xxx_4);
xxx_5  = xxx==5;
xxx_5  = double(xxx_5);
mmm_1  = mmm==1;
mmm_1  = double(mmm_1);
mmm_2  = mmm==2;
mmm_2  = double(mmm_2);
mmm_3  = mmm==3;
mmm_3  = double(mmm_3);
mmm_4  = mmm==4;
mmm_4  = double(mmm_4);
mmm_5  = mmm==5;
mmm_5  = double(mmm_5);
mmm_6  = mmm==6;
mmm_6  = double(mmm_6);
mmm_7  = mmm==7;
mmm_7  = double(mmm_7);
mmm_8  = mmm==8;
mmm_8  = double(mmm_8);
mmm_9  = mmm==9;
mmm_9  = double(mmm_9);
mmm_10 = mmm==10;
mmm_10 = double(mmm_10);
mmm_11 = mmm==11;
mmm_11 = double(mmm_11);
mmm_12 = mmm==12;
mmm_12 = double(mmm_12);
mmm_13 = mmm==13;
mmm_13 = double(mmm_13);
mmm_14 = mmm==14;
mmm_14 = double(mmm_14);
mmm_15 = mmm==15;
mmm_15 = double(mmm_15);
mmm_16 = mmm==16;
mmm_16 = double(mmm_16);
mmm_17 = mmm==17;
mmm_17 = double(mmm_17);
mmm_18 = mmm==18;
mmm_18 = double(mmm_18);
mmm_19 = mmm==19;
mmm_19 = double(mmm_19);
mmm_20 = mmm==20;
mmm_20 = double(mmm_20);
ttt_1  = ttt==1;
ttt_1  = double(ttt_1);
ttt_2  = ttt==2;
ttt_2  = double(ttt_2);
ttt_3  = ttt==3;
ttt_3  = double(ttt_3);
ttt_4  = ttt==4;
ttt_4  = double(ttt_4);
ttt_5  = ttt==5;
ttt_5  = double(ttt_5);
ttt_6  = ttt==6;
ttt_6  = double(ttt_6);
ttt_7  = ttt==7;
ttt_7  = double(ttt_7);
ttt_8  = ttt==8;
ttt_8  = double(ttt_8);
ttt_9  = ttt==9;
ttt_9  = double(ttt_9);
ttt_10  = ttt==10;
ttt_10  = double(ttt_10);

% I left j =5 ; m=20 ; t=10 outside to get away from multicollinearity:

xxxx    = [xxx_1 xxx_2 xxx_3 xxx_4 xxx_5];
mmmm    = [mmm_2 mmm_3 mmm_4 mmm_5 mmm_6 mmm_7 mmm_8 mmm_9 mmm_10 mmm_11 mmm_12 mmm_13 mmm_14 mmm_15 mmm_16 mmm_17 mmm_18 mmm_19 mmm_20];
tttt    = [ttt_2 ttt_3 ttt_4 ttt_5 ttt_6 ttt_7 ttt_8 ttt_9 ttt_10];

X_xi    = [xxxx mmmm tttt];                                                % Stacking regressors together

end