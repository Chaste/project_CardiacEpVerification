
mono1 = [0.1, 0.04729256677, 0.2559326838;
0.05, 0.01219747854, 0.1241794945;
0.025, 0.003076150583, 0.06157125564;
0.0125, 0.0007707704265, 0.030719;
0.00625, 0.0001928018545, 0.01535111097;];

mono2 = [0.1, 0.09725475116, 0.9603208893;
0.05, 0.02701293006, 0.4708969787;
0.025, 0.006950391271, 0.2337494882;
0.0125, 0.001750637086, 0.1166425258;
0.00625, 0.0004385267801, 0.05829111062;];

mono3 = [];

bi1 = [0.1, 0.04265058497, 0.2537848512, 0.02988107835, 0.1723212284;
0.05, 0.01089005107, 0.1238652902, 0.00760416821, 0.08669301031;
0.025, 0.002739314674, 0.06153026531, 0.00191092449, 0.04339679202;
0.0125, 0.0006859231324, 0.03071381977, 0.0004783746599, 0.02170398909;
0.00625, 0.000171550361, 0.01535046176, 0.0001196347713, 0.01085267045;];

bi2 = [0.1, 0.06757124444, 0.9378114481, 0.04767141994, 0.6344914118;
0.05, 0.01775058515, 0.4670640583, 0.0125103767, 0.3266563191;
0.025, 0.004496896232, 0.2332270697, 0.00316834263, 0.1644645911;
0.0125, 0.001128029025, 0.1165725359, 0.0007946991365, 0.0823727219;
0.00625, 0.0002822467237, 0.0582810478, 0.0001988392088, 0.04120385967;];

bi3 = [];

bath2 = [0.1, 0.04404732625, 0.2364742997, 0.05985167155, 0.1633438949;
0.05, 0.01134407088, 0.1137769974, 0.01533843483, 0.08013854566;
0.025, 0.002860772769, 0.05627417462, 0.003864296406, 0.03978798753;
0.0125, 0.0007168090283, 0.02805656314, 0.0009681450983, 0.01984697711;
0.00625, 0.0001793041513, 0.01401773733, 0.0002421824826, 0.009915136989;];


do_mono = 0;
do_bi_V = 0; 
do_bi_phi_e = 0;
do_bath_V = 0;
do_bath_phi_e = 1;

do_bath = (do_bath_V || do_bath_phi_e);


if(do_mono + do_bi_V + do_bi_phi_e + do_bath_V + do_bath_phi_e ~= 1)
  error('Bad options');
end;

%% monodomain
if(do_mono)
  dat1 = mono1;
  dat2 = mono2;
  dat3 = mono3;
end;

%% bidomain voltage
if(do_bi_V)
  dat1 = bi1(:,1:3);
  dat2 = bi2(:,1:3);
  dat3 = bi3(:,1:3);
end;

%% bidomain phi_e
if(do_bi_phi_e)
  dat1 = bi1(:,[1,4:5]);
  dat2 = bi2(:,[1,4:5]);
  dat3 = bi3(:,[1,4:5]);
end;

if(~do_bath)

  for ind=2:3
    if(ind==2)
       type = 'k*-';
    else
       type = 'k*--';
    end
    loglog(dat1(:,1), dat1(:,ind), type);
    if(ind==2)
      hold on;
      %pos = get(gca,'pos'); pos([3 4]) = pos([3 4])*0.5;set(gca,'pos',pos);
    end;

    loglog(dat2(:,1), dat2(:,ind), type);
    loglog(dat3(:,1), dat3(:,ind), type);

    if(ind==2)
      text(dat1(1,1)/0.9,dat1(1,ind),'1D');
      text(dat2(1,1)/0.9,dat2(1,ind),'2D');
      text(dat3(1,1)/0.9,dat3(1,ind),'3D');
    end;
  end

  ind = 3;
  x1 = dat1(2,1); x2 = dat1(4,1);
  loglog([x1 x2], 100*[x1 x2], 'k');
  loglog([x2 x2 x1], 100*[x2 x1 x1], 'k');
  
  ind = 2;
  x1 = dat1(2,1); x2 = dat1(4,1);
  y1 = dat1(2,ind)*(dat1(2,ind)/dat2(2,ind)); y2 = dat1(3,ind)*(dat1(3,ind)/dat2(3,ind));
  loglog([x1 x2], [x1 x2].^2, 'k');
  loglog([x2 x1 x1], [x2 x2 x1].^2, 'k');

  axis([10^(-2.5) 10^-0.5 1e-4 10^(1)]);
  xlabel('Spatial stepsize, h');

  if(do_mono)
    ylabel('Error');
  end;
  if(do_bi_V)
    ylabel('Error in voltage');
  end;
  if(do_bi_phi_e)
    ylabel('Error in extracellular potential');
  end;
end;


if(do_bath)
  if(do_bath_V) 
    dat = bath2(:,[1,2:3]);
  else
    dat = bath2(:,[1,4:5]);
  end

  figure
  loglog(dat(:,1), dat(:,2), 'k*-');
  hold on;
  %pos = get(gca,'pos'); pos([3 4]) = pos([3 4])*0.5;set(gca,'pos',pos);
  loglog(dat(:,1), dat(:,3), 'k*--');

  x1 = dat(2,1); x2 = dat(4,1);
  loglog([x1 x2], 10*[x1 x2], 'k');
  loglog([x2 x2 x1], 10*[x2 x1 x1], 'k');

  loglog([x1 x2], [x1 x2].^2, 'k');
  loglog([x2 x1 x1], [x2 x2 x1].^2, 'k');

  axis([10^(-2.5) 10^-0.5 1e-4 10^0]);

  xlabel('Spatial stepsize, h');
  if(do_bath_V)
    ylabel('Error in voltage');
  else
    ylabel('Error in extracellular potential');
  end;
end;


