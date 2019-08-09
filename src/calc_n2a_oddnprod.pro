pro calc_n2a_oddnprod,zmajnow,zei,chema,chemq,p_ave,p_n2d,p_no

P7_2D  =zei.PEI_n2*zmajnow.n2den*.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.zz-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.oden+chemq.q_81(3)*zmajnow.o2den+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.zz-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.zz-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.zz-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]];+$

P7_2D=p7_2d+zei.aEI_n2*zmajnow.n2den*.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.zz-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.o+chemq.q_81(3)*zmajnow.o2+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.zz-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.zz-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.zz-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]];+$

p_n2d=p7_2d

; 3,4,5,6 are vibrational levels of N2A which are fit independantly
; p_ave=fit parameters (constant from gaussian scaling of N2A production per N2+ production)  
;********WITHOUT TEMPERATURE DEPENDENCE****************
P4_NO = zei.PEI_n2*zmajnow.n2den*0.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.z-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.oden+chemq.q_81(3)*zmajnow.o2den+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.z-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.z-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.z-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]]

P4_NO = P4_NO+zei.aEI_n2*zmajnow.n2den*0.57*zmajnow.oden*[[[[1./[[p_ave(0,3)+p_ave(1,3)*exp((-1.*(zmajnow.z-p_ave(2,3))^2)/(2.*p_ave(3,3)^2))]]]/[chemq.q_80(3)*zmajnow.oden+chemq.q_81(3)*zmajnow.o2den+chema.A_8(3)]]*chemq.q_80(3)]+$
[[[1./[[p_ave(0,4)+p_ave(1,4)*exp((-1.*(zmajnow.z-p_ave(2,4))^2)/(2.*p_ave(3,4)^2))]]]/[chemq.q_80(4)*zmajnow.oden+chemq.q_81(4)*zmajnow.o2den+chema.A_8(4)]]*chemq.q_80(4)]+$
[[[1./[[p_ave(0,5)+p_ave(1,5)*exp((-1.*(zmajnow.z-p_ave(2,5))^2)/(2.*p_ave(3,5)^2))]]]/[chemq.q_80(5)*zmajnow.oden+chemq.q_81(5)*zmajnow.o2den+chema.A_8(5)]]*chemq.q_80(5)]+$
[[[1./[[p_ave(0,6)+p_ave(1,6)*exp((-1.*(zmajnow.z-p_ave(2,6))^2)/(2.*p_ave(3,6)^2))]]]/[chemq.q_80(6)*zmajnow.oden+chemq.q_81(6)*zmajnow.o2den+chema.A_8(6)]]*chemq.q_80(6)]]

p_no=p4_no

return
end