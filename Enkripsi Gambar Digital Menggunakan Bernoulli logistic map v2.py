from Tkinter import *
from tkFileDialog import *
from PIL import Image
from PIL import ImageTk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
import math
import os
import timeit

def gantiframe(halaman):
	halaman.tkraise()

####################
# chaos

def sebarBoL(r, l, x0):
	nilai = r*l*x0*(1 - x0)
	lantai = np.floor(nilai)
	return nilai - lantai
def bifur(rku, minl, maksl):
	
	sebarku = sebarBoL

	a.clear()
	a.set_title('Diagram bifurkasi dengan r = %s' %rku)
	n = 10000
	lamda = np.linspace(minl, maksl, n)

	iterasi = 1000
	last = 100

	x0ku = 1e-5 * np.ones(n)

	for i in range(iterasi):
		# x0 = Bern(lamda,x0)
		# x0 = Log(lamda, x0)
		x0ku = sebarku(rku, lamda, x0ku)
		# x0 = LogBern(10, lamda, x0)

		if i >= (iterasi - last):
			a.plot(lamda, x0ku, 'b,', alpha = 0.04)
	a.set_xlabel('lambda')
	a.set_ylabel('x')
	f.tight_layout()
	canvas.show()

def kunciBoL(n, urutan, r, l, x0):
	kunci = np.zeros(n, dtype=np.uint8)
	
	for i in range(urutan):
		nilai = r*l*x0*(1 - x0)
		lantai = np.floor(nilai)
		x0 = nilai - lantai

	for i in range(n):
		nilai = r*l*x0*(1 - x0)
		lantai = np.floor(nilai)
		x0 = nilai - lantai
		kunci[i] = x0 * 1000 % 256
	return kunci

def plot(nku, urku, rku, lku, x0ku):
	a.clear()
	a.set_title('Sebaran %s nilai x pertama' %nku)
	a.set_xlabel('n')
	a.set_ylabel('x_n')
	a.plot(kunciBoL(nku, urku, rku, lku, x0ku),'b',
		kunciBoL(nku, urku, rku, lku, x0ku),'k.')
	# a.plot(kunciBoL(banyakku, mulaiku, lku, rku, x0ku),'k.')
	f.tight_layout()
	canvas.show()


# cobweb

def bol(urutan, r, l, x0, banyak):
	res = np.zeros(banyak)
	x = x0
	for i in range(urutan):
		x = r*l*x*(1-x)-np.floor(r*l*x*(1-x))
	for i in range(banyak):
		res[i] = x
		x = r*l*x*(1-x)-np.floor(r*l*x*(1-x))
	return res

def cobweb3(r,l, f, x0, banyak):
	a.clear()
	maxi = int(r*l/4.)+1
	x = [0,1]
	for i in range(maxi):
		if r*l == 0:
			tambahan = 0.5
		else:
			tambahan = np.sqrt(1./4 - i/float(r*l))
		print tambahan
		x.append(0.5 - tambahan)
		x.append(0.5 + tambahan)
	x.sort()
	print x
	for i in range(len(x)-1):
		cvec = np.arange(x[i]+0.00001,x[i+1], 0.00001)
		a.plot(cvec, r*l*cvec*(1-cvec)-np.floor(r*l*cvec*(1-cvec)),'b')
	cvec = np.arange(0,1,0.01)
	a.plot(cvec, cvec, linestyle='--')
	x = bol(0, r, l, x0, banyak)
	for i in range(1, banyak):
		a.plot([x[i-1],x[i-1]],[x[i-1],x[i]], color='gray')
		a.plot([x[i-1],x[i]],[x[i],x[i]], color='gray')
	a.set_title('cobweb')
	a.set_xlabel('x_n')
	a.set_ylabel('x_(n+1)')
	f.tight_layout()
	canvas.show()



########################
# nist
def bil2bin(vektor):
	b = [format(i,'08b') for i in vektor]
	return ''.join(b)
def tulisfile(namafile, barisanbin):
	fileku = open(namafile, 'w')
	fileku.write(barisanbin)
	fileku.close()
def kunci2txt(nku, urku, rku, lku, x0ku):
	namafile = asksaveasfilename()
	bil = kunciBoL(nku, urku, rku, lku, x0ku)
	binku = bil2bin(bil)
	tulisfile(namafile, binku)

# monobit
def monobitku(databin):
	vektor = [ int(i) for i in databin]
	n = len(vektor)

	x = [2*i -1 for i in vektor]
	absx = abs(sum(x))
	sobs = absx / np.sqrt(n)
	pvalue = math.erfc(sobs/np.sqrt(2))
	if pvalue > 0.01:
		keputusan = 'acak (monobit)'
	else:
		keputusan = 'takacak (monobit)'
	return n, absx, sobs, pvalue, keputusan

# runtest
def runtest(databin):
	vektor = [int(i) for i in databin]
	n = len(vektor)
	pi = 1.0 * sum(vektor)/n
	tau = 2/np.sqrt(n)
	abspimintau = abs(pi-0.5)
	vobs = len(databin.replace('0', ' ').split()) + len(databin.replace('1' , ' ').split())
	if abspimintau >= tau:
		pvalue = 0
	else:
		pvalue = math.erfc(abs(vobs-2*n*pi*(1-pi)) / (2 * pi * (1 - pi) * np.sqrt(2*n)))
	if pvalue > 0.01:
		keputusan = 'acak (runtest)'
	else:
		keputusan = 'takacak (runtest)'
	return n, pi, tau, abspimintau, vobs, pvalue, keputusan

def kunci2tampil(nku, urku, rku, lku, x0ku):
	bil = kunciBoL(nku, urku, rku, lku, x0ku)
	binku = bil2bin(bil)
	n, absx, sobs, pvaluemono, kepmono = monobitku(binku)
	nmonon.set(n)
	nmonoabsx.set(absx)
	nmonosobs.set(sobs)
	nmonopvalue.set(pvaluemono)
	nmonokeputusan.set(kepmono)

# def kunci2run(nku, urku, rku, lku, x0ku):
# 	bil = kunciBoL(nku, urku, rku, lku, x0ku)
# 	binku = bil2bin(bil)
	n, pi, tau, abspimintau, vobs, pvaluerun, keprun = runtest(binku)
	nrunn.set(n)
	nrunpi.set(pi)
	nruntau.set(tau)
	nrunabspimintau.set(abspimintau)
	nrunvobs.set(vobs)
	nrunpvalue.set(pvaluerun)
	nrunkeputusan.set(keprun)




###############################
# kripsi
############################

def simpangb(almentry, alm):
	formatku = [('portable network graphics','*.png'),('semua','*')]
	gb = Image.open(alm)
	namafile = asksaveasfilename(parent=root, filetypes=formatku, title='simpan gambar sebagai...')
	gb.save(namafile)
	almentry.set(namafile)

def resizeku(gb):
	pasli, lasli = gb.size
	pa = float(pasli)
	la = float(lasli)

	p2 = 350
	l2 = p2 / pa * la
	if l2 > 250:
		l2 = 250
		p2 = l2 / la * pa
	p2 = int(p2)
	l2 = int(l2)
	return gb.resize((p2,l2))

# proses
def kripsi(gb, ur, r, l, x0, mode):
	mgb = np.asarray(gb)
	if len(mgb.shape) == 2: 
		d = 1			# kalau filenya grayscale atau format .gif
		b,k = mgb.shape
	else:
		b,k,d = mgb.shape

	vgb = mgb.reshape(b*k*d)

	kunciku = kunciBoL(b*k*d, ur, r, l, x0)
	vcipher = np.zeros(b*k*d,dtype = np.uint8)
	vcipher[0] = vgb[0] ^ kunciku[0]
	if mode == 'en':
		for i in range(1,b*k*d):
			vcipher[i] = vgb[i] ^ kunciku[i] ^ vcipher[i-1]
	else:
		for i in range(1,b*k*d):
			vcipher[i] = vgb[i] ^ kunciku[i] ^ vgb[i-1]
	if d == 1:			# kalau filenya grayscale atau format .gif
		mcipher = vcipher.reshape(b,k)
	else:
		mcipher = vcipher.reshape(b,k,d)
	# gbcipher = Image.fromarray(mcipher)
	# gbcipher.show()
	# gbcipher.save('cip2.png')
	return Image.fromarray(mcipher)

# -------------
# enkripsi
# -----------------
def CariP():
	formatku = [('jpeg','*.jpg'), ('portable network graphics','*.png'),('semua','*')]
	alamat = askopenfilename(parent=root, filetypes=formatku)
	ealm0.set(alamat)
	ukuran = str(float(os.path.getsize(alamat)) / 1000) + 'kb'
	# eket0 ukuranPlainteks.set('ukuran: ' + ukuran)
	gb = Image.open(alamat)
	pasli, lasli = gb.size
	eket0.set(ukuran + ', ' + str(pasli) + 'x' + str(lasli) + 'px')
	
	gbresize = resizeku(gb)
	grTk = ImageTk.PhotoImage(gbresize)
	egb0.configure(image=grTk)
	egb0.image = grTk

def enk(alm, ur, r, l, x0):
	awalwaktu = timeit.default_timer()
	gambar = Image.open(alm)
	gbcipku = kripsi(gambar, ur, r, l, x0, 'en')
	waktu = timeit.default_timer() - awalwaktu
	ewaktu.set('%.5f detik' %waktu )

	gbcipku.save('cipku.png')
	pcip,lcip = gbcipku.size
	ukuran = str(float(os.path.getsize('cipku.png')) / 1000)+ 'kb'
	ekett.set(ukuran + ', ' + str(pcip) + 'x' + str(lcip) + 'px')
	
	res = resizeku(gbcipku)
	gbcipku2 = ImageTk.PhotoImage(res)
	egbt.configure(image=gbcipku2)
	egbt.image = gbcipku2
	ealmt.set('cipku.png')

# def esimpan():
# 	gb = Image.open('cipku.png')
# 	namafile = asksaveasfilename()
# 	gb.save(namafile)





###############################
# dekripsi
############################
def CariC():
	formatku = [('portable network graphics','*.png'),('semua','*')]
	alamat = askopenfilename(parent=root, filetypes=formatku)
	dalm0.set(alamat)
	ukuran = str(float(os.path.getsize(alamat)) / 1000) + 'kb'
	# eket0 ukuranPlainteks.set('ukuran: ' + ukuran)
	gb = Image.open(alamat)
	pasli, lasli = gb.size
	dket0.set(ukuran + ', ' + str(pasli) + 'x' + str(lasli) + 'px')
	
	gbresize = resizeku(gb)
	grTk = ImageTk.PhotoImage(gbresize)
	dgb0.configure(image=grTk)
	dgb0.image = grTk

def dek(alm, ur, r, l, x0):
	awalwaktu = timeit.default_timer()
	gb0 = Image.open(alm)
	gbt = kripsi(gb0, ur, r, l, x0, 'de')
	waktu = timeit.default_timer() - awalwaktu
	dwaktu.set('%.5f detik' %waktu )

	gbt.save('plainku.png')
	pcip,lcip = gbt.size
	ukuran = str(float(os.path.getsize('plainku.png')) / 1000)+ 'kb'
	dkett.set(ukuran + ', ' + str(pcip) + 'x' + str(lcip) + 'px')
	
	res = resizeku(gbt)
	gbt2 = ImageTk.PhotoImage(res)
	dgbt.configure(image=gbt2)
	dgbt.image = gbt2
	dalmt.set('plainku.png')

# def dsimpan():
# 	gb = Image.open('plainku.png')
# 	namafile = asksaveasfilename()
# 	gb.save(namafile)


######################
# histogram
# -----------------
def CariH(almentry, gbku):
	formatku = [('semua','*'), ('jpeg','*.jpg'), ('portable network graphics','*.png')]
	alamat = askopenfilename(parent=root, filetypes=formatku)
	almentry.set(alamat)
	
	gb = Image.open(alamat)
	pasli, lasli = gb.size
	gbresize = resizeku(gb)
	grTk = ImageTk.PhotoImage(gbresize)
	gbku.configure(image=grTk)
	gbku.image = grTk


def histoku(almku, hgbhis, cvku):
	gb = Image.open(almku)
	mgb = np.asarray(gb)
	if len(mgb.shape) == 2:
		grayhis = gb.histogram()
		hgbhis.clear()
		hispgray = hgbhis.add_subplot(111)
		hispgray.bar(range(256),grayhis, edgecolor = "none")
		hispgray.set_title('Histogram Grayscale')
		hispgray.set_xlim([0,255])
	else:
		r,g,b = mgb[:,:,0], mgb[:,:,1], mgb[:,:,2]
		rku = Image.fromarray(r)
		gku = Image.fromarray(g)
		bku = Image.fromarray(b)
		rhis = rku.histogram()
		ghis = gku.histogram()
		bhis = bku.histogram()
		# plot 
		hgbhis.clear()
		hispr = hgbhis.add_subplot(311)
		hispr.bar(range(256),rhis, edgecolor = "none")
		hispr.set_title('Histogram r')
		hispr.set_xlim([0,255])
		
		hispg = hgbhis.add_subplot(312)
		# plt.subplot(3,1,2)
		hispg.bar(range(256),ghis, edgecolor = "none")
		hispg.set_title('Histogram g')
		hispg.set_xlim([0,255])
		
		hispb = hgbhis.add_subplot(313)
		# plt.subplot(3,1,3)
		hispb.bar(range(256),bhis, edgecolor = "none")
		hispb.set_title('Histogram b')
		hispb.set_xlim([0,255])
		# plt.savefig('histogram rgb.png')
		# plt.show()
	hgbhis.tight_layout()
	cvku.show()

def hproses():
	histoku(halm0.get(), hgbhis0, hcv0)
	histoku(halmt.get(), hgbhist, hcvt)


#################
# psnr
# --------------
def psnr(almgbA,almgbB):
	ga = Image.open(almgbA)
	gb = Image.open(almgbB)
	ma = np.asarray(ga)
	mb = np.asarray(gb)
	b1,k1,d1 = ma.shape
	b2,k2,d2 = mb.shape

# print b1,k1,d1
# print b2,k2,d2

	if b1==b2 and k1==k2 and d1==d2:
		n = b1*k1*d1
		mra = ma.reshape(b1*k1*d1)
		mrb = mb.reshape(b1*k1*d1)
		mse = 1./n*sum((mra - mrb)**2)
		print mse
		if mse == 0:
			psnr = 'tak terdefinisi'
		else:
			psnr = 20*np.log10(255/np.sqrt(mse))
		print psnr
	else:
		mse = 'tak dapat dihitung, ukuran gambar beda'
		psnr = 'tak dapat dihitung, ukuran gambar beda'
		print 'ukuran tidak sama'
	return mse, psnr

def tampilPSNR():
	mseku, psnrku = psnr(psalm0.get(), psalmt.get())
	psmse.set(mseku)
	pspsnr.set(psnrku)











##############################


# ini untuk frame dasar, tempat semua frame
root = Tk()
root.title('mebl')
dasar = Frame(root)
dasar.pack(expand=True)
dasar.grid_rowconfigure(0, weight=1)
dasar.grid_columnconfigure(0, weight=1)


fAwal = Frame(dasar)
fEnk = Frame(dasar)
fDek = Frame(dasar)
fAn = Frame(dasar)
fHis = Frame(dasar)
fNIST = Frame(dasar)
fPSNR = Frame(dasar)
fPar = Frame(dasar)

for F in(fAwal, fEnk, fDek, fAn, fHis, fNIST, fPSNR, fPar):
	F.grid(row=0, column=0, sticky='news')
###########################

# ini untuk menu
menuku = Menu(root)
root.config(menu=menuku)

pilihMenu = Menu(menuku)
menuku.add_cascade(label='Menu', menu=pilihMenu)
pilihMenu.add_command(label='Enkripsi', command=lambda:gantiframe(fEnk))
pilihMenu.add_command(label='Dekripsi', command=lambda:gantiframe(fDek))
pilihMenu.add_command(label='Halaman Utama', command=lambda:gantiframe(fAwal))

parMenu = Menu(menuku)
menuku.add_cascade(label='Parameter', menu=parMenu)
parMenu.add_command(label='Diag.Bifurkasi, Plot sebaran, & Cobweb', command=lambda:gantiframe(fPar))
parMenu.add_command(label='Uji NIST', command=lambda:gantiframe(fNIST))

anMenu = Menu(menuku)
menuku.add_cascade(label='Analisis', menu=anMenu)
anMenu.add_command(label='Histogram', command=lambda:gantiframe(fHis))
anMenu.add_command(label='Uji NIST', command=lambda:gantiframe(fNIST))
anMenu.add_command(label='PSNR', command=lambda:gantiframe(fPSNR))
########################


# ini di frame Awal
Label(fAwal, text='Enkripsi Gambar Digital dengan\nBernoulli Logistic Map').pack()
render = Image.open('gblatex/bol.png')
bolawal = ImageTk.PhotoImage(render)
Label(fAwal, image=bolawal).pack()
Label(fAwal, text='Indra Bayu Muktyas').pack()



# ini di frame fPar parameter
# cv = Canvas(fPar, width=300, height=200, bg='yellow')
# cv.pack(expand=True, fill='both')

# gbku = Image.open('gblatex/cobweblagi2.png')
# gifku = ImageTk.PhotoImage(gbku)
# cv.create_image(50,10, image=gifku, anchor='nw')

#####################
# parameter global
vurutan = IntVar()
vr = DoubleVar()
vl = DoubleVar()
vx0 = DoubleVar()
#######################

par = Frame(fPar)
par.pack(fill='both', expand=True)


Label(par, text='urutan:').grid(row=1, column=0, sticky='e')
Label(par, text='r:').grid(row=2, column=0, sticky='e')
Label(par, text='lambda:').grid(row=3, column=0, sticky='e')
Label(par, text='x0:').grid(row=4, column=0, sticky='e')
urutan = Entry(par, textvariable=vurutan)
r = Entry(par, textvariable=vr)
l = Entry(par, textvariable=vl)
x0 = Entry(par, textvariable=vx0)
urutan.grid(row=1, column=1)
r.grid(row=2, column=1)
l.grid(row=3, column=1)
x0.grid(row=4, column=1)



Label(par, text='min lambda:').grid(row=1, column=3, sticky='e')
Label(par, text='maks lambda:').grid(row=2, column=3, sticky='e')

minlambda=DoubleVar()
makslambda=DoubleVar()
n=IntVar()
minlambda.set(0)
makslambda.set(5)
n.set(200)
Entry(par, textvariable=minlambda).grid(row=1, column=4)
Entry(par, textvariable=makslambda).grid(row=2, column=4)

# fplot = LabelFrame(par)
# fplot.grid(row=3, column=2, rowspan=2, columnspan=3, sticky='news')
Label(par, text='n:').grid(row=4, column=3, sticky='e')
Entry(par, textvariable=n).grid(row=4, column=4)
# tombol diag.bifur dan plot
Button(par, text='diag.bifur', command=lambda:bifur(vr.get(), 
	minlambda.get(), makslambda.get())).grid(row=1, column=2, rowspan=2)
Button(par, text='plot x_n', command=lambda:plot(n.get(),
	vurutan.get(), vr.get(), vl.get(), vx0.get())).grid(row=3, column=2, rowspan=2, sticky='es')
Button(par, text='cobweb', command=lambda:cobweb3(
	vr.get(), vl.get(), f, vx0.get(),n.get())).grid(row=1, column=5, rowspan=2)
# cobweb3(r,l, ax, x0, banyak=20

f = Figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)
# t = np.arange(0.0,3.0,0.01)
# s = np.sin(2*np.pi*t)

# a.plot(t,s)


# a tk.DrawingArea
canvas = FigureCanvasTkAgg(f, master=fPar)
canvas.show()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2TkAgg( canvas, fPar )
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)



##################################
# ini frame enkripsi
###############################

parEnk = Frame(fEnk)
parEnk.pack(fill='both', expand=True)


Label(parEnk, text='urutan:').grid(row=1, column=0, sticky='e')
Label(parEnk, text='r:').grid(row=2, column=0, sticky='e')
Label(parEnk, text='lambda:').grid(row=3, column=0, sticky='e')
Label(parEnk, text='x0:').grid(row=4, column=0, sticky='e')
urutan = Entry(parEnk, textvariable=vurutan)
r = Entry(parEnk, textvariable=vr)
l = Entry(parEnk, textvariable=vl)
x0 = Entry(parEnk, textvariable=vx0)
urutan.grid(row=1, column=1)
r.grid(row=2, column=1)
l.grid(row=3, column=1)
x0.grid(row=4, column=1)
Button(parEnk, text='parameter\nlanjutan', 
	command=lambda:gantiframe(fPar)).grid(row=1, column=2, rowspan=4)

prosesEnk = Frame(fEnk)
prosesEnk.pack(fill='both', expand=True)

# plainteks
ealm0 = StringVar()
Label(prosesEnk, text='Plainteks').grid(row=0, column=0)
Entry(prosesEnk, textvariable=ealm0).grid(row=1, column=0)
ecari = Button(prosesEnk, text='cari', command=CariP).grid(row=1, column=1)

eket0 = StringVar()
egb0 = Label(prosesEnk, text='Belum ada gambar', relief='sunken')
egb0.grid(row=2, column=0, rowspan=3, columnspan=2, sticky='news')
Label(prosesEnk, text='ket.', textvariable=eket0).grid(row=5, column=0, columnspan=2)

ewaktu = DoubleVar()
etblEnk = Button(prosesEnk, text='enkripsi>>', command=lambda:enk(
	ealm0.get(), vurutan.get(), vr.get(), vl.get(), vx0.get())).grid(row=3, column=2)

Label(prosesEnk, text='waktu :').grid(row=4, column=2)
Entry(prosesEnk, textvariable=ewaktu).grid(row=5, column=2)

# cipherteks
ealmt = StringVar()
Label(prosesEnk, text='Ciperteks').grid(row=0, column=3)
Entry(prosesEnk, textvariable=ealmt).grid(row=1, column=3, columnspan=2, sticky='we')

ekett = StringVar()
egbt = Label(prosesEnk, text='Belum ada gambar', relief='sunken')
egbt.grid(row=2, column=3, rowspan=3, columnspan=2, sticky='news')
Label(prosesEnk, text='ket.', textvariable=ekett).grid(row=5, column=3, columnspan=2)
etblSimpan = Button(prosesEnk, text='simpan', 
	command=lambda:simpangb(ealmt, 'cipku.png')).grid(row=6, column=4, sticky='es')





# ini frame dekripsi

parDek = Frame(fDek)
parDek.pack(fill='both', expand=True)


Label(parDek, text='urutan:').grid(row=1, column=0, sticky='e')
Label(parDek, text='r:').grid(row=2, column=0, sticky='e')
Label(parDek, text='lambda:').grid(row=3, column=0, sticky='e')
Label(parDek, text='x0:').grid(row=4, column=0, sticky='e')
urutan = Entry(parDek, textvariable=vurutan)
r = Entry(parDek, textvariable=vr)
l = Entry(parDek, textvariable=vl)
x0 = Entry(parDek, textvariable=vx0)
urutan.grid(row=1, column=1)
r.grid(row=2, column=1)
l.grid(row=3, column=1)
x0.grid(row=4, column=1)
Button(parDek, text='parameter\nlanjutan', 
	command=lambda:gantiframe(fPar)).grid(row=1, column=2, rowspan=4)

prosesDek = Frame(fDek)
prosesDek.pack(fill='both', expand=True)

# cipherteks
dalm0 = StringVar()
Label(prosesDek, text='Cipherteks').grid(row=0, column=0)
Entry(prosesDek, textvariable=dalm0).grid(row=1, column=0)
dcari = Button(prosesDek, text='cari', command=CariC).grid(row=1, column=1)

dket0 = StringVar()
dgb0 = Label(prosesDek, text='Belum ada gambar', relief='sunken')
dgb0.grid(row=2, column=0, rowspan=3, columnspan=2, sticky='news')
Label(prosesDek, text='ket.', textvariable=dket0).grid(row=5, column=0, columnspan=2)

dwaktu = DoubleVar()
dtblDe = Button(prosesDek, text='dekripsi>>', command=lambda:dek(
	dalm0.get(), vurutan.get(), vr.get(), vl.get(), vx0.get())).grid(row=3, column=2)
Label(prosesDek, text='waktu :').grid(row=4, column=2)
Entry(prosesDek, textvariable=dwaktu).grid(row=5, column=2)

# decipherteks
dalmt = StringVar()
Label(prosesDek, text='Gambar terdekripsi').grid(row=0, column=3)
Entry(prosesDek, textvariable=dalmt).grid(row=1, column=3, columnspan=2, sticky='we')

dkett = StringVar()
dgbt = Label(prosesDek, text='Belum ada gambar', relief='sunken')
dgbt.grid(row=2, column=3, rowspan=3, columnspan=2, sticky='news')
Label(prosesDek, text='ket.', textvariable=dkett).grid(row=5, column=3, columnspan=2)
dtblSimpan = Button(prosesDek, text='simpan', 
	command=lambda:simpangb(dalmt, 'plainku.png')).grid(row=6, column=4, sticky='es')




##############################
# histogram
halm0 = StringVar()
Label(fHis, text='Gambar awal').grid(row=0, column=0)
Entry(fHis, textvariable=halm0).grid(row=1, column=0)
hcari0 = Button(fHis, text='cari', command=lambda:CariH(
	halm0, hgb0)).grid(row=1, column=1)
hgb0 = Label(fHis, text='belum ada gambar', relief='sunken')
hgb0.grid(row=2, column=0, columnspan=2, sticky='news')

halmt = StringVar()
Label(fHis, text='Gambar akhir').grid(row=0, column=3)
Entry(fHis, textvariable=halmt).grid(row=1, column=3)
hcarit = Button(fHis, text='cari', command=lambda:CariH(
	halmt, hgbt)).grid(row=1, column=4)
hgbt = Label(fHis, text='belum ada gambar', relief='sunken')
hgbt.grid(row=2, column=3, columnspan=2, sticky='news')

htbl = Button(fHis, text='proses\nvv',command=hproses).grid(row=2, column=2)


# tampilan histogramnya
hgbhis0= Figure(figsize=(4,4), dpi=100)
# hisp = gbHisPlainteks.add_subplot(311)
hcv0 = FigureCanvasTkAgg(hgbhis0, fHis)
hcv0.get_tk_widget().grid(row=3, column=0, columnspan=2)

# hisCipherteks = Button(frameAnalisis, text='histogram', command=histokucip)
# hisCipherteks.grid(row=1, column=3, sticky='w')
# harusnya pakai canvas buat nggambar histogramnya
hgbhist = Figure(figsize=(4,4), dpi=100)
# hisc = gbHisCipherteks.add_subplot(311)
hcvt = FigureCanvasTkAgg(hgbhist, fHis)
hcvt.get_tk_widget().grid(row=3, column=3, columnspan=2)


#######################
# PSNR

psalm0 = StringVar()
Label(fPSNR, text='Gambar awal').grid(row=0, column=0)
Entry(fPSNR, textvariable=psalm0).grid(row=1, column=0)

psgb0 = Label(fPSNR, text='belum ada gambar', relief='sunken')
psgb0.grid(row=2, column=0, columnspan=2, sticky='news')
pscari0 = Button(fPSNR, text='cari', command=lambda:CariH(psalm0,
	psgb0)).grid(row=1, column=1)

psalmt = StringVar()
Label(fPSNR, text='Gambar terdekripsi').grid(row=0, column=3)
Entry(fPSNR, textvariable=psalmt).grid(row=1, column=3)

psgbt = Label(fPSNR, text='belum ada gambar', relief='sunken')
psgbt.grid(row=2, column=3, columnspan=2, sticky='news')
pscarit = Button(fPSNR, text='cari', command=lambda:CariH(psalmt,
	psgbt)).grid(row=1, column=4)
pstbl = Button(fPSNR, text='proses\nvv', command=tampilPSNR).grid(row=2, column=2)

psmse = DoubleVar()
pspsnr = DoubleVar()
psjadi = DoubleVar()
Label(fPSNR, text='MSE =').grid(row=3, column=0)
Label(fPSNR, text='PSNR =').grid(row=4, column=0)
Label(fPSNR, textvariable=psjadi)

Entry(fPSNR, textvariable=psmse).grid(row=3, column=1)
Entry(fPSNR, textvariable=pspsnr).grid(row=4, column=1)



########################
# uji NIST
nbuat = LabelFrame(fNIST, text='Buat bitstring txt')
nbuat.grid(row=1, column=0)

Label(nbuat, text='urutan:').grid(row=1, column=0, sticky='e')
Label(nbuat, text='r:').grid(row=2, column=0, sticky='e')
Label(nbuat, text='lambda:').grid(row=3, column=0, sticky='e')
Label(nbuat, text='x0:').grid(row=4, column=0, sticky='e')
Label(nbuat, text='banyaknya x:').grid(row=5, column=0, sticky='e')

urutan = Entry(nbuat, textvariable=vurutan)
r = Entry(nbuat, textvariable=vr)
l = Entry(nbuat, textvariable=vl)
x0 = Entry(nbuat, textvariable=vx0)
banyak = Entry(nbuat, textvariable=n)

urutan.grid(row=1, column=1)
r.grid(row=2, column=1)
l.grid(row=3, column=1)
x0.grid(row=4, column=1)
banyak.grid(row=5, column=1)
nsimpanbit = Button(nbuat, text='simpan', command=lambda:kunci2txt(n.get(),
	vurutan.get(), vr.get(), vl.get(), vx0.get())).grid(row=6, column=1, sticky='e')
nproses = Button(fNIST, text='proses', command=lambda:kunci2tampil(
	n.get(), vurutan.get(), vr.get(), vl.get(), vx0.get() )).grid(row=2, column=0)


# frequency (monobit) test
nfrek = LabelFrame(fNIST, text='frequency (monobit) test')
nfrek.grid(row=3, column=0)
Label(nfrek, text='n =').grid(row=1, column=0, sticky='e')
Label(nfrek, text='absx =').grid(row=2, column=0, sticky='e')
Label(nfrek, text='sobs =').grid(row=3, column=0, sticky='e')
Label(nfrek, text='pvalue =').grid(row=4, column=0, sticky='e')
Label(nfrek, text='keputusan =').grid(row=5, column=0, sticky='e')

nmonon = IntVar()
nmonoabsx = IntVar()
nmonosobs = DoubleVar()
nmonopvalue = DoubleVar()
nmonokeputusan = StringVar()
Entry(nfrek, textvariable=nmonon).grid(row=1, column=1, sticky='w')
Entry(nfrek, textvariable=nmonoabsx).grid(row=2, column=1, sticky='w')
Entry(nfrek, textvariable=nmonosobs).grid(row=3, column=1, sticky='w')
Entry(nfrek, textvariable=nmonopvalue).grid(row=4, column=1, sticky='w')
Entry(nfrek, textvariable=nmonokeputusan).grid(row=5, column=1, sticky='w')
# -----------

# run test
nrun = LabelFrame(fNIST, text='run test')
nrun.grid(row=4, column=0)
Label(nrun, text='n =').grid(row=1, column=0, sticky='e')
Label(nrun, text='pi =').grid(row=2, column=0, sticky='e')
Label(nrun, text='tau =').grid(row=3, column=0, sticky='e')
Label(nrun, text='|pi-tau| =').grid(row=4, column=0, sticky='e')
Label(nrun, text='vobs =').grid(row=5, column=0, sticky='e')
Label(nrun, text='pvalue =').grid(row=6, column=0, sticky='e')
Label(nrun, text='keputusan =').grid(row=7, column=0, sticky='e')

nrunn = IntVar()
nrunpi = DoubleVar()
nruntau = DoubleVar()
nrunabspimintau = DoubleVar()
nrunvobs = DoubleVar()
nrunpvalue = DoubleVar()
nrunkeputusan = StringVar()
Entry(nrun, textvariable=nrunn).grid(row=1, column=1, sticky='w')
Entry(nrun, textvariable=nrunpi).grid(row=2, column=1, sticky='w')
Entry(nrun, textvariable=nruntau).grid(row=3, column=1, sticky='w')
Entry(nrun, textvariable=nrunabspimintau).grid(row=4, column=1, sticky='w')
Entry(nrun, textvariable=nrunvobs).grid(row=5, column=1, sticky='w')
Entry(nrun, textvariable=nrunpvalue).grid(row=6, column=1, sticky='w')
Entry(nrun, textvariable=nrunkeputusan).grid(row=7, column=1, sticky='w')
# -----------


gantiframe(fAwal)
root.mainloop()