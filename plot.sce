subplot(331)
xtitle("u")
xgrid()
plot(hoge.time,hoge.values(:,1))
subplot(332)
xtitle("v")
xgrid()
plot(hoge.time,hoge.values(:,2))
subplot(333)
xtitle("w")
xgrid()
plot(hoge.time,hoge.values(:,3))

subplot(334)
xtitle("p")
xgrid()
plot(hoge.time,hoge.values(:,4))
subplot(335)
xtitle("q")
xgrid()
plot(hoge.time,hoge.values(:,5))
subplot(336)
xtitle("r")
xgrid()
plot(hoge.time,hoge.values(:,6))

subplot(337)
xtitle("phi")
xgrid()
plot(hoge.time,hoge.values(:,7))
subplot(338)
xtitle("theta")
xgrid()
plot(hoge.time,hoge.values(:,8))
subplot(339)
xtitle("psi")
xgrid()
plot(hoge.time,hoge.values(:,9))

//Input
scf()
subplot(411)
xtitle("d1")
xgrid()
plot(delta.time,delta.values(:,1))
subplot(412)
xtitle("d2")
xgrid()
plot(delta.time,delta.values(:,2))
subplot(413)
xtitle("d3")
xgrid()
plot(delta.time,delta.values(:,3))
subplot(414)
xtitle("d4")
xgrid()
plot(delta.time,delta.values(:,4))
