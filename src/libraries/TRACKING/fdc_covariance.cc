

// The following was auto-generated from the gen_covariance_code.C
// macro using a ROOT file generated with the fdc_covariance_tree
// plugin.


double GetFDCCovariance(int layer1, int layer2);
double GetFDCCathodeCovariance(int layer1, int layer2);


//-------------------------
// GetFDCCovariance
//-------------------------
double GetFDCCovariance(int layer1, int layer2)
{
	if(layer1<1 || layer2>24 || layer2<1 || layer2>24)return 0.0;
	if(layer2<layer1){
		int tmp = layer1;
		layer1 = layer2;
		layer2 = tmp;
	}

	switch(layer1){
		case 1:
			if(layer2==1)return 0.000585634;
			if(layer2==2)return 0.000585637;
			if(layer2==3)return 0.000585505;
			if(layer2==4)return 0.000580329;
			if(layer2==5)return 0.000577621;
			if(layer2==6)return 0.000577471;
			if(layer2==7)return 0.000337295;
			if(layer2==8)return 0.000265735;
			if(layer2==9)return 0.000281526;
			if(layer2==10)return 0.000220917;
			if(layer2==11)return 0.000199922;
			if(layer2==12)return 9.48379e-05;
			if(layer2==13)return -0.000621933;
			if(layer2==14)return -0.00085377;
			if(layer2==15)return -0.00100895;
			if(layer2==16)return -0.00101167;
			if(layer2==17)return -0.00109397;
			if(layer2==18)return -0.00122267;
			if(layer2==19)return -0.002534;
			if(layer2==20)return -0.0026614;
			if(layer2==21)return -0.00297017;
			if(layer2==22)return -0.00294615;
			if(layer2==23)return -0.00307489;
			if(layer2==24)return -0.00327564;
			break; // layer 1
		case 2:
			if(layer2==2)return 0.000588198;
			if(layer2==3)return 0.000586648;
			if(layer2==4)return 0.000583233;
			if(layer2==5)return 0.000571174;
			if(layer2==6)return 0.0005714;
			if(layer2==7)return 0.000323606;
			if(layer2==8)return 0.000241671;
			if(layer2==9)return 0.000241261;
			if(layer2==10)return 0.000209708;
			if(layer2==11)return 0.000207312;
			if(layer2==12)return 0.000104978;
			if(layer2==13)return -0.000697916;
			if(layer2==14)return -0.000897349;
			if(layer2==15)return -0.00102077;
			if(layer2==16)return -0.0010817;
			if(layer2==17)return -0.00110823;
			if(layer2==18)return -0.00126598;
			if(layer2==19)return -0.00259498;
			if(layer2==20)return -0.00264318;
			if(layer2==21)return -0.00302039;
			if(layer2==22)return -0.00292141;
			if(layer2==23)return -0.00320805;
			if(layer2==24)return -0.00328668;
			break; // layer 2
		case 3:
			if(layer2==3)return 0.000610334;
			if(layer2==4)return 0.000576509;
			if(layer2==5)return 0.000580319;
			if(layer2==6)return 0.000552877;
			if(layer2==7)return 0.000270971;
			if(layer2==8)return 0.000269639;
			if(layer2==9)return 0.000313385;
			if(layer2==10)return 0.000204131;
			if(layer2==11)return 0.000213183;
			if(layer2==12)return 0.000142326;
			if(layer2==13)return -0.0006886;
			if(layer2==14)return -0.000836541;
			if(layer2==15)return -0.000867524;
			if(layer2==16)return -0.00089191;
			if(layer2==17)return -0.00118841;
			if(layer2==18)return -0.00114409;
			if(layer2==19)return -0.00262366;
			if(layer2==20)return -0.00276578;
			if(layer2==21)return -0.00302496;
			if(layer2==22)return -0.00274118;
			if(layer2==23)return -0.00329716;
			if(layer2==24)return -0.00322495;
			break; // layer 3
		case 4:
			if(layer2==4)return 0.000663642;
			if(layer2==5)return 0.00057822;
			if(layer2==6)return 0.000603196;
			if(layer2==7)return 0.000508455;
			if(layer2==8)return 0.000383883;
			if(layer2==9)return 0.000420884;
			if(layer2==10)return 0.000469223;
			if(layer2==11)return 0.00029155;
			if(layer2==12)return 0.00028835;
			if(layer2==13)return -0.000305435;
			if(layer2==14)return -0.000420861;
			if(layer2==15)return -0.000793265;
			if(layer2==16)return -0.000887853;
			if(layer2==17)return -0.00069152;
			if(layer2==18)return -0.000970917;
			if(layer2==19)return -0.00212543;
			if(layer2==20)return -0.00233911;
			if(layer2==21)return -0.00263463;
			if(layer2==22)return -0.00286518;
			if(layer2==23)return -0.00275972;
			if(layer2==24)return -0.00311307;
			break; // layer 4
		case 5:
			if(layer2==5)return 0.000762076;
			if(layer2==6)return 0.000578621;
			if(layer2==7)return 0.000344426;
			if(layer2==8)return 0.000501357;
			if(layer2==9)return 0.000475929;
			if(layer2==10)return 0.000359146;
			if(layer2==11)return 0.000263931;
			if(layer2==12)return 0.000297406;
			if(layer2==13)return -0.000500437;
			if(layer2==14)return -0.000267775;
			if(layer2==15)return -0.000595586;
			if(layer2==16)return -0.000903103;
			if(layer2==17)return -0.000563898;
			if(layer2==18)return -0.000759776;
			if(layer2==19)return -0.00246909;
			if(layer2==20)return -0.00263376;
			if(layer2==21)return -0.00224664;
			if(layer2==22)return -0.00282813;
			if(layer2==23)return -0.00312209;
			if(layer2==24)return -0.00318908;
			break; // layer 5
		case 6:
			if(layer2==6)return 0.00091888;
			if(layer2==7)return 0.000592297;
			if(layer2==8)return 0.000361344;
			if(layer2==9)return 0.000600295;
			if(layer2==10)return 0.000481805;
			if(layer2==11)return 0.000418612;
			if(layer2==12)return 0.000475262;
			if(layer2==13)return -0.000160005;
			if(layer2==14)return -0.00033444;
			if(layer2==15)return -0.000574285;
			if(layer2==16)return -0.000972138;
			if(layer2==17)return -0.000565014;
			if(layer2==18)return -0.000808768;
			if(layer2==19)return -0.00186498;
			if(layer2==20)return -0.00204238;
			if(layer2==21)return -0.00240424;
			if(layer2==22)return -0.0025992;
			if(layer2==23)return -0.00244965;
			if(layer2==24)return -0.00274029;
			break; // layer 6
		case 7:
			if(layer2==7)return 0.00941613;
			if(layer2==8)return 0.00243384;
			if(layer2==9)return 0.00197798;
			if(layer2==10)return 0.00198357;
			if(layer2==11)return 0.00159895;
			if(layer2==12)return 0.00217037;
			if(layer2==13)return 0.00540402;
			if(layer2==14)return 0.00346997;
			if(layer2==15)return 0.00354863;
			if(layer2==16)return 0.00232144;
			if(layer2==17)return 0.00755751;
			if(layer2==18)return 0.0039261;
			if(layer2==19)return 0.00384286;
			if(layer2==20)return 0.00944355;
			if(layer2==21)return 0.00565062;
			if(layer2==22)return 0.00136788;
			if(layer2==23)return 0.00990618;
			if(layer2==24)return 0.00620554;
			break; // layer 7
		case 8:
			if(layer2==8)return 0.012477;
			if(layer2==9)return 0.00336062;
			if(layer2==10)return 0.00236514;
			if(layer2==11)return 0.00236415;
			if(layer2==12)return 0.00322261;
			if(layer2==13)return 0.00475278;
			if(layer2==14)return 0.00672238;
			if(layer2==15)return 0.0101226;
			if(layer2==16)return 0.00384073;
			if(layer2==17)return 0.0133649;
			if(layer2==18)return 0.0110647;
			if(layer2==19)return 0.000943711;
			if(layer2==20)return 0.0119423;
			if(layer2==21)return 0.0172422;
			if(layer2==22)return 0.00180998;
			if(layer2==23)return 0.00971831;
			if(layer2==24)return 0.014592;
			break; // layer 8
		case 9:
			if(layer2==9)return 0.0109299;
			if(layer2==10)return 0.00263205;
			if(layer2==11)return 0.00191202;
			if(layer2==12)return 0.00214378;
			if(layer2==13)return 0.00508611;
			if(layer2==14)return 0.0037997;
			if(layer2==15)return 0.00686033;
			if(layer2==16)return 0.00407503;
			if(layer2==17)return 0.00651406;
			if(layer2==18)return 0.00784606;
			if(layer2==19)return 0.00432215;
			if(layer2==20)return 0.00561573;
			if(layer2==21)return 0.00929072;
			if(layer2==22)return 0.00430716;
			if(layer2==23)return 0.00506619;
			if(layer2==24)return 0.00863102;
			break; // layer 9
		case 10:
			if(layer2==10)return 0.0126246;
			if(layer2==11)return 0.00205085;
			if(layer2==12)return 0.00266086;
			if(layer2==13)return 0.00694037;
			if(layer2==14)return 0.00696621;
			if(layer2==15)return 0.00301847;
			if(layer2==16)return 0.00527808;
			if(layer2==17)return 0.00657487;
			if(layer2==18)return 0.00332406;
			if(layer2==19)return 0.00717002;
			if(layer2==20)return 0.00689147;
			if(layer2==21)return 0.00488396;
			if(layer2==22)return 0.00461135;
			if(layer2==23)return 0.00663182;
			if(layer2==24)return 0.0028906;
			break; // layer 10
		case 11:
			if(layer2==11)return 0.0137119;
			if(layer2==12)return 0.00202584;
			if(layer2==13)return 0.00307288;
			if(layer2==14)return 0.00566122;
			if(layer2==15)return 0.00451167;
			if(layer2==16)return 0.00243922;
			if(layer2==17)return 0.00643825;
			if(layer2==18)return 0.00528795;
			if(layer2==19)return 0.00257297;
			if(layer2==20)return 0.00759791;
			if(layer2==21)return 0.00873083;
			if(layer2==22)return 0.00224811;
			if(layer2==23)return 0.00626871;
			if(layer2==24)return 0.00829066;
			break; // layer 11
		case 12:
			if(layer2==12)return 0.0149375;
			if(layer2==13)return 0.00669672;
			if(layer2==14)return 0.00532517;
			if(layer2==15)return 0.00923196;
			if(layer2==16)return 0.00514177;
			if(layer2==17)return 0.00858525;
			if(layer2==18)return 0.00916853;
			if(layer2==19)return 0.00660234;
			if(layer2==20)return 0.00686408;
			if(layer2==21)return 0.0117095;
			if(layer2==22)return 0.00690335;
			if(layer2==23)return 0.00704723;
			if(layer2==24)return 0.0109506;
			break; // layer 12
		case 13:
			if(layer2==13)return 0.0424831;
			if(layer2==14)return 0.0132986;
			if(layer2==15)return 0.0112832;
			if(layer2==16)return 0.0114952;
			if(layer2==17)return 0.0191486;
			if(layer2==18)return 0.014654;
			if(layer2==19)return 0.0265102;
			if(layer2==20)return 0.0251338;
			if(layer2==21)return 0.0204195;
			if(layer2==22)return 0.0229189;
			if(layer2==23)return 0.0266061;
			if(layer2==24)return 0.019977;
			break; // layer 13
		case 14:
			if(layer2==14)return 0.0487339;
			if(layer2==15)return 0.0129352;
			if(layer2==16)return 0.0113439;
			if(layer2==17)return 0.0210415;
			if(layer2==18)return 0.0162345;
			if(layer2==19)return 0.0162787;
			if(layer2==20)return 0.0296003;
			if(layer2==21)return 0.029972;
			if(layer2==22)return 0.0144352;
			if(layer2==23)return 0.0306181;
			if(layer2==24)return 0.026856;
			break; // layer 14
		case 15:
			if(layer2==15)return 0.0547796;
			if(layer2==16)return 0.0134318;
			if(layer2==17)return 0.0266401;
			if(layer2==18)return 0.0293692;
			if(layer2==19)return 0.0171676;
			if(layer2==20)return 0.0272509;
			if(layer2==21)return 0.0500961;
			if(layer2==22)return 0.023925;
			if(layer2==23)return 0.0278054;
			if(layer2==24)return 0.0489375;
			break; // layer 15
		case 16:
			if(layer2==16)return 0.0481183;
			if(layer2==17)return 0.0145537;
			if(layer2==18)return 0.015112;
			if(layer2==19)return 0.0260147;
			if(layer2==20)return 0.0241747;
			if(layer2==21)return 0.0201543;
			if(layer2==22)return 0.0321924;
			if(layer2==23)return 0.0241489;
			if(layer2==24)return 0.0211445;
			break; // layer 16
		case 17:
			if(layer2==17)return 0.0724264;
			if(layer2==18)return 0.0326953;
			if(layer2==19)return 0.0166213;
			if(layer2==20)return 0.0551036;
			if(layer2==21)return 0.0552371;
			if(layer2==22)return 0.0154214;
			if(layer2==23)return 0.0549375;
			if(layer2==24)return 0.0523486;
			break; // layer 17
		case 18:
			if(layer2==18)return 0.0680758;
			if(layer2==19)return 0.0219335;
			if(layer2==20)return 0.0321661;
			if(layer2==21)return 0.0626003;
			if(layer2==22)return 0.0276625;
			if(layer2==23)return 0.0324177;
			if(layer2==24)return 0.0611823;
			break; // layer 18
		case 19:
			if(layer2==19)return 0.103584;
			if(layer2==20)return 0.0381741;
			if(layer2==21)return 0.0361168;
			if(layer2==22)return 0.0609626;
			if(layer2==23)return 0.0433971;
			if(layer2==24)return 0.0380554;
			break; // layer 19
		case 20:
			if(layer2==20)return 0.127328;
			if(layer2==21)return 0.0655345;
			if(layer2==22)return 0.0363031;
			if(layer2==23)return 0.0968891;
			if(layer2==24)return 0.0702118;
			break; // layer 20
		case 21:
			if(layer2==21)return 0.146928;
			if(layer2==22)return 0.042056;
			if(layer2==23)return 0.0665459;
			if(layer2==24)return 0.114058;
			break; // layer 21
		case 22:
			if(layer2==22)return 0.115602;
			if(layer2==23)return 0.0406957;
			if(layer2==24)return 0.0476336;
			break; // layer 22
		case 23:
			if(layer2==23)return 0.146012;
			if(layer2==24)return 0.0737714;
			break; // layer 23
		case 24:
			if(layer2==24)return 0.158888;
			break; // layer 24
	} // switch for layer1 25

	return 0.0;
}


// NOTE: At this point, the following is incorrect. It does not
// Properly account for the Lorentz deflections that are in the
// simulated data file. It's better to use the above for the
// MULS error along the wire as well for the time being.
//    April 20, 2009  DL


//-------------------------
// GetFDCCathodeCovariance
//-------------------------
double GetFDCCathodeCovariance(int layer1, int layer2)
{
	if(layer1<1 || layer2>24 || layer2<1 || layer2>24)return 0.0;
	if(layer2<layer1){
		int tmp = layer1;
		layer1 = layer2;
		layer2 = tmp;
	}

	switch(layer1){
		case 1:
			if(layer2==1)return 0.000259212;
			if(layer2==2)return -1.48748e-06;
			if(layer2==3)return 1.49637e-06;
			if(layer2==4)return -3.95846e-06;
			if(layer2==5)return -2.62431e-06;
			if(layer2==6)return -7.89397e-07;
			if(layer2==7)return 2.85425e-05;
			if(layer2==8)return 1.10816e-05;
			if(layer2==9)return 1.69028e-05;
			if(layer2==10)return 1.82918e-05;
			if(layer2==11)return -1.4731e-05;
			if(layer2==12)return 1.48368e-05;
			if(layer2==13)return 5.5183e-05;
			if(layer2==14)return 1.10052e-06;
			if(layer2==15)return 3.33326e-05;
			if(layer2==16)return 4.88196e-05;
			if(layer2==17)return 4.25319e-05;
			if(layer2==18)return 4.02736e-05;
			if(layer2==19)return 8.80668e-05;
			if(layer2==20)return 0.00011175;
			if(layer2==21)return 8.40024e-05;
			if(layer2==22)return 8.43182e-05;
			if(layer2==23)return 5.43089e-05;
			if(layer2==24)return 4.01294e-05;
			break; // layer 1
		case 2:
			if(layer2==2)return 0.000257954;
			if(layer2==3)return -2.17538e-07;
			if(layer2==4)return 3.05871e-06;
			if(layer2==5)return 1.05325e-05;
			if(layer2==6)return -2.02414e-05;
			if(layer2==7)return 1.69923e-05;
			if(layer2==8)return 8.42202e-05;
			if(layer2==9)return -0.00010199;
			if(layer2==10)return -8.60192e-06;
			if(layer2==11)return 9.5151e-05;
			if(layer2==12)return -0.000119313;
			if(layer2==13)return -5.02435e-05;
			if(layer2==14)return 0.000146529;
			if(layer2==15)return -0.000198696;
			if(layer2==16)return -7.37034e-05;
			if(layer2==17)return 0.000119301;
			if(layer2==18)return -0.000206375;
			if(layer2==19)return -0.000194452;
			if(layer2==20)return 0.000123303;
			if(layer2==21)return -0.000278509;
			if(layer2==22)return -0.000156683;
			if(layer2==23)return 0.000115377;
			if(layer2==24)return -0.000255137;
			break; // layer 2
		case 3:
			if(layer2==3)return 0.000282503;
			if(layer2==4)return 2.49987e-05;
			if(layer2==5)return -3.31924e-05;
			if(layer2==6)return 7.65871e-05;
			if(layer2==7)return 0.00021057;
			if(layer2==8)return -0.000127246;
			if(layer2==9)return 0.000350758;
			if(layer2==10)return 0.000240732;
			if(layer2==11)return -0.000122853;
			if(layer2==12)return 0.000408539;
			if(layer2==13)return 0.000508711;
			if(layer2==14)return -2.98293e-05;
			if(layer2==15)return 0.00046555;
			if(layer2==16)return 0.000503079;
			if(layer2==17)return 6.94061e-05;
			if(layer2==18)return 0.000487515;
			if(layer2==19)return 0.000742624;
			if(layer2==20)return 0.000241953;
			if(layer2==21)return 0.000493597;
			if(layer2==22)return 0.000790171;
			if(layer2==23)return 0.000201507;
			if(layer2==24)return 0.000524748;
			break; // layer 3
		case 4:
			if(layer2==4)return 0.000344896;
			if(layer2==5)return 6.70214e-05;
			if(layer2==6)return 9.35823e-05;
			if(layer2==7)return 0.000658416;
			if(layer2==8)return 0.000560364;
			if(layer2==9)return 0.000199577;
			if(layer2==10)return 0.000817453;
			if(layer2==11)return 0.000727325;
			if(layer2==12)return 0.000190464;
			if(layer2==13)return 0.00123537;
			if(layer2==14)return 0.00138458;
			if(layer2==15)return 4.15432e-05;
			if(layer2==16)return 0.00126318;
			if(layer2==17)return 0.00151974;
			if(layer2==18)return -5.78469e-05;
			if(layer2==19)return 0.00099527;
			if(layer2==20)return 0.00202511;
			if(layer2==21)return -0.000791871;
			if(layer2==22)return 0.00139361;
			if(layer2==23)return 0.00180254;
			if(layer2==24)return -0.00117376;
			break; // layer 4
		case 5:
			if(layer2==5)return 0.000480153;
			if(layer2==6)return -0.000176782;
			if(layer2==7)return 0.000348708;
			if(layer2==8)return 0.00136624;
			if(layer2==9)return -0.00114986;
			if(layer2==10)return 0.000233608;
			if(layer2==11)return 0.00152214;
			if(layer2==12)return -0.00137914;
			if(layer2==13)return -0.00040315;
			if(layer2==14)return 0.00213068;
			if(layer2==15)return -0.00211505;
			if(layer2==16)return -0.000439276;
			if(layer2==17)return 0.00201734;
			if(layer2==18)return -0.002279;
			if(layer2==19)return -0.00176081;
			if(layer2==20)return 0.00219088;
			if(layer2==21)return -0.00349883;
			if(layer2==22)return -0.00169445;
			if(layer2==23)return 0.00186074;
			if(layer2==24)return -0.00390196;
			break; // layer 5
		case 6:
			if(layer2==6)return 0.000700908;
			if(layer2==7)return 0.00111843;
			if(layer2==8)return -0.000725816;
			if(layer2==9)return 0.00217654;
			if(layer2==10)return 0.00166731;
			if(layer2==11)return -0.000618979;
			if(layer2==12)return 0.00249643;
			if(layer2==13)return 0.00354854;
			if(layer2==14)return 8.93556e-05;
			if(layer2==15)return 0.00324484;
			if(layer2==16)return 0.00372352;
			if(layer2==17)return 0.000449007;
			if(layer2==18)return 0.00333059;
			if(layer2==19)return 0.00526369;
			if(layer2==20)return 0.00147875;
			if(layer2==21)return 0.00358997;
			if(layer2==22)return 0.00537658;
			if(layer2==23)return 0.00146612;
			if(layer2==24)return 0.00306215;
			break; // layer 6
		case 7:
			if(layer2==7)return 0.0119196;
			if(layer2==8)return 0.00537461;
			if(layer2==9)return 0.0074455;
			if(layer2==10)return 0.0120928;
			if(layer2==11)return 0.00762934;
			if(layer2==12)return 0.00816088;
			if(layer2==13)return 0.0251839;
			if(layer2==14)return 0.0171282;
			if(layer2==15)return 0.00711942;
			if(layer2==16)return 0.0264409;
			if(layer2==17)return 0.021432;
			if(layer2==18)return 0.00620119;
			if(layer2==19)return 0.0293231;
			if(layer2==20)return 0.0329347;
			if(layer2==21)return -0.00264357;
			if(layer2==22)return 0.0309255;
			if(layer2==23)return 0.0327569;
			if(layer2==24)return -0.0059781;
			break; // layer 7
		case 8:
			if(layer2==8)return 0.0127684;
			if(layer2==9)return -0.00804617;
			if(layer2==10)return 0.00582225;
			if(layer2==11)return 0.0152876;
			if(layer2==12)return -0.0100867;
			if(layer2==13)return 0.00246289;
			if(layer2==14)return 0.023448;
			if(layer2==15)return -0.0191363;
			if(layer2==16)return 0.00186881;
			if(layer2==17)return 0.0236892;
			if(layer2==18)return -0.0212029;
			if(layer2==19)return -0.00898381;
			if(layer2==20)return 0.0254384;
			if(layer2==21)return -0.0331911;
			if(layer2==22)return -0.00820424;
			if(layer2==23)return 0.0231141;
			if(layer2==24)return -0.037125;
			break; // layer 8
		case 9:
			if(layer2==9)return 0.0176829;
			if(layer2==10)return 0.00818921;
			if(layer2==11)return -0.007947;
			if(layer2==12)return 0.0210759;
			if(layer2==13)return 0.0275057;
			if(layer2==14)return -0.00597435;
			if(layer2==15)return 0.0307075;
			if(layer2==16)return 0.0299251;
			if(layer2==17)return -0.00121394;
			if(layer2==18)return 0.0325854;
			if(layer2==19)return 0.0454171;
			if(layer2==20)return 0.00923176;
			if(layer2==21)return 0.0356939;
			if(layer2==22)return 0.0467138;
			if(layer2==23)return 0.0117531;
			if(layer2==24)return 0.0341138;
			break; // layer 9
		case 10:
			if(layer2==10)return 0.0159691;
			if(layer2==11)return 0.00824885;
			if(layer2==12)return 0.00907031;
			if(layer2==13)return 0.0282879;
			if(layer2==14)return 0.0206978;
			if(layer2==15)return 0.00801328;
			if(layer2==16)return 0.0296676;
			if(layer2==17)return 0.0235688;
			if(layer2==18)return 0.00721606;
			if(layer2==19)return 0.0325955;
			if(layer2==20)return 0.0348687;
			if(layer2==21)return -0.0021366;
			if(layer2==22)return 0.033933;
			if(layer2==23)return 0.0341879;
			if(layer2==24)return -0.00765239;
			break; // layer 10
		case 11:
			if(layer2==11)return 0.018272;
			if(layer2==12)return -0.0100808;
			if(layer2==13)return 0.00710489;
			if(layer2==14)return 0.0300333;
			if(layer2==15)return -0.0215207;
			if(layer2==16)return 0.00631188;
			if(layer2==17)return 0.030934;
			if(layer2==18)return -0.0238843;
			if(layer2==19)return -0.00585993;
			if(layer2==20)return 0.0350924;
			if(layer2==21)return -0.0394446;
			if(layer2==22)return -0.0045144;
			if(layer2==23)return 0.0321007;
			if(layer2==24)return -0.0448891;
			break; // layer 11
		case 12:
			if(layer2==12)return 0.0251355;
			if(layer2==13)return 0.0326002;
			if(layer2==14)return -0.00883425;
			if(layer2==15)return 0.0376816;
			if(layer2==16)return 0.0355568;
			if(layer2==17)return -0.00350749;
			if(layer2==18)return 0.0403662;
			if(layer2==19)return 0.0562171;
			if(layer2==20)return 0.00928686;
			if(layer2==21)return 0.0455012;
			if(layer2==22)return 0.0573382;
			if(layer2==23)return 0.013248;
			if(layer2==24)return 0.0434772;
			break; // layer 12
		case 13:
			if(layer2==13)return 0.0746805;
			if(layer2==14)return 0.0314059;
			if(layer2==15)return 0.0435293;
			if(layer2==16)return 0.0837659;
			if(layer2==17)return 0.0438281;
			if(layer2==18)return 0.0430851;
			if(layer2==19)return 0.111447;
			if(layer2==20)return 0.0816691;
			if(layer2==21)return 0.0275505;
			if(layer2==22)return 0.113626;
			if(layer2==23)return 0.0861752;
			if(layer2==24)return 0.0179491;
			break; // layer 13
		case 14:
			if(layer2==14)return 0.0634629;
			if(layer2==15)return -0.0300999;
			if(layer2==16)return 0.0309582;
			if(layer2==17)return 0.0687258;
			if(layer2==18)return -0.0358981;
			if(layer2==19)return 0.0150537;
			if(layer2==20)return 0.0882759;
			if(layer2==21)return -0.0726655;
			if(layer2==22)return 0.0162644;
			if(layer2==23)return 0.0862753;
			if(layer2==24)return -0.0853443;
			break; // layer 14
		case 15:
			if(layer2==15)return 0.0740832;
			if(layer2==16)return 0.0505133;
			if(layer2==17)return -0.0240409;
			if(layer2==18)return 0.0806006;
			if(layer2==19)return 0.0899995;
			if(layer2==20)return -0.00857435;
			if(layer2==21)return 0.0961208;
			if(layer2==22)return 0.0974814;
			if(layer2==23)return 0.000700254;
			if(layer2==24)return 0.0988446;
			break; // layer 15
		case 16:
			if(layer2==16)return 0.0892268;
			if(layer2==17)return 0.0442551;
			if(layer2==18)return 0.050842;
			if(layer2==19)return 0.123553;
			if(layer2==20)return 0.0862985;
			if(layer2==21)return 0.0351894;
			if(layer2==22)return 0.12626;
			if(layer2==23)return 0.0926381;
			if(layer2==24)return 0.0252649;
			break; // layer 16
		case 17:
			if(layer2==17)return 0.0786674;
			if(layer2==18)return -0.0306203;
			if(layer2==19)return 0.033929;
			if(layer2==20)return 0.107187;
			if(layer2==21)return -0.0739739;
			if(layer2==22)return 0.0334797;
			if(layer2==23)return 0.107798;
			if(layer2==24)return -0.0866908;
			break; // layer 17
		case 18:
			if(layer2==18)return 0.0893095;
			if(layer2==19)return 0.0971627;
			if(layer2==20)return -0.0154624;
			if(layer2==21)return 0.113637;
			if(layer2==22)return 0.103139;
			if(layer2==23)return -0.0082286;
			if(layer2==24)return 0.114272;
			break; // layer 18
		case 19:
			if(layer2==19)return 0.191898;
			if(layer2==20)return 0.0949705;
			if(layer2==21)return 0.0985812;
			if(layer2==22)return 0.202249;
			if(layer2==23)return 0.112417;
			if(layer2==24)return 0.0882125;
			break; // layer 19
		case 20:
			if(layer2==20)return 0.173644;
			if(layer2==21)return -0.0796099;
			if(layer2==22)return 0.0958044;
			if(layer2==23)return 0.182096;
			if(layer2==24)return -0.0976008;
			break; // layer 20
		case 21:
			if(layer2==21)return 0.182626;
			if(layer2==22)return 0.104415;
			if(layer2==23)return -0.0722736;
			if(layer2==24)return 0.190527;
			break; // layer 21
		case 22:
			if(layer2==22)return 0.211122;
			if(layer2==23)return 0.110645;
			if(layer2==24)return 0.0941412;
			break; // layer 22
		case 23:
			if(layer2==23)return 0.194591;
			if(layer2==24)return -0.091654;
			break; // layer 23
		case 24:
			if(layer2==24)return 0.206504;
			break; // layer 24
	} // switch for layer1 25

	return 0.0;
}
