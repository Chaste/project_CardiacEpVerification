% note, tested in octave not matlab


function draw_convergence_plots

handle = figure;
hold on;

subplot(3,2,1);
do_one_plot(1);

for i=2:5
  subplot(3,2,i+1);
  do_one_plot(i);
end;

FN = findall(handle,'-property','FontName');
set(FN,'FontName','Arial')

FN = findall(handle,'-property','LineWidth');
set(FN,'LineWidth',1)




function do_one_plot(index)

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

mono3 = [0.1, 0.1850259458, 2.083578;
0.05, 0.05883684687, 1.055421835;
0.025, 0.01571919181, 0.5268933745;
0.0125, 0.004102424571, 0.2632805022;];

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

bi3 = [0.1, 0.1167727055, 1.99370425, 0.08253275914, 1.350677342;
0.05, 0.03257310981, 1.037230142, 0.02301724628, 0.7258621174;
0.025, 0.008372181866, 0.5243032404, 0.005915685911, 0.3697861914;
0.0125, 0.00210780985, 0.2629039327, 0.001489328919, 0.1857819454;];

bath2 = [0.1, 0.04404732625, 0.2364742997, 0.05985167155, 0.1633438949;
0.05, 0.01134407088, 0.1137769974, 0.01533843483, 0.08013854566;
0.025, 0.002860772769, 0.05627417462, 0.003864296406, 0.03978798753;
0.0125, 0.0007168090283, 0.02805656314, 0.0009681450983, 0.01984697711;
0.00625, 0.0001793041513, 0.01401773733, 0.0002421824826, 0.009915136989;];


do_mono = (index==1);
do_bi_V = (index==2); 
do_bi_phi_e = (index==3);
do_bath_V = (index==4);
do_bath_phi_e = (index==5);

do_bath = (do_bath_V || do_bath_phi_e);


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
    loglog(dat1(:,1), dat1(:,ind), type, "markersize", 3);
    if(ind==2)
      hold on;
    end;

    loglog(dat2(:,1), dat2(:,ind), type, "markersize", 3);
    loglog(dat3(:,1), dat3(:,ind), type, "markersize", 3);

    if(ind==2)
      text(dat1(1,1)*0.7,dat1(1,ind)/4,'1D','fontsize',8);
      text(dat2(5,1)/0.9,dat2(5,ind)*3,'2D','fontsize',8);
      text(dat3(4,1)/0.9,dat3(4,ind)*3,'3D','fontsize',8);
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
  if(index>3)
    xlabel('Spatial stepsize, h');
  end;

  if(do_mono)
    title('M-1D, M-2D, M-3D');
    ylabel('Error');
  end;
  if(do_bi_V)
    title('B-1D, B-2D, B-3D');
    ylabel('Error in V');
  end;
  if(do_bi_phi_e)
    title('B-1D, B-2D, B-3D');
    ylabel('Error in \phi_e');
  end;
end;


if(do_bath)
  if(do_bath_V) 
    dat = bath2(:,[1,2:3]);
  else
    dat = bath2(:,[1,4:5]);
  end

  loglog(dat(:,1), dat(:,2), 'k*-', "markersize", 3);
  hold on;
  loglog(dat(:,1), dat(:,3), 'k*--', "markersize", 3);

  x1 = dat(2,1); x2 = dat(4,1);
  loglog([x1 x2], 10*[x1 x2], 'k');
  loglog([x2 x2 x1], 10*[x2 x1 x1], 'k');

  loglog([x1 x2], [x1 x2].^2, 'k');
  loglog([x2 x1 x1], [x2 x2 x1].^2, 'k');

  axis([10^(-2.5) 10^-0.5 1e-4 10^0]);

  xlabel('Spatial stepsize, h');
  if(do_bath_V)
    title('BB-2D*');
    ylabel('Error in V');
  else
    title('BB-2D*');
    ylabel('Error in \phi_e');
  end;
end;


