using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace TestImageAnalog
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            Random r = new Random();
            double fe = 4000000;
            double SNR = 30;
            int echantParPixel = 3;
            double[,] depart = ImageLecture((Bitmap)Bitmap.FromFile("D:/lab/imtestaller.bmp"), echantParPixel,1);
            Console.WriteLine("T transmission : "+depart.GetLength(1)/fe+" s");
            Console.WriteLine("T par pixel : " + echantParPixel / fe + " s");
            Console.WriteLine("F transmission : " + fe/depart.GetLength(1) + " Hz");
            double[] sig = Transmettre(depart, fe);
            sig = BruiterSignal(sig, SNR, fe, ref r);
            //enregistrerMatrice(vectVersMat(DSP(sig,sig.Length)), "D:/lab/imspectre.csv");
            double[,] arrivee = Recevoir(sig,5,100,fe);
            ImageEcriture(arrivee,170,100,0.2).Save("D:/lab/imtestretour.bmp");
            System.Diagnostics.Process.Start("D:/lab/imtestretour.bmp");
        }
        public static double[,] ImageLecture(Bitmap image, int nbRepet, int res)
        {
            int n = 0;
            for (int x = 0; x < image.Width; x += res)
            {
                for (int y = 0; y < image.Height; y += res)
                {
                    n++;
                }
            }
            double[,] retour = new double[5, n * nbRepet];
            double r ;
            double v ;
            double b ;
            double px ;
            double py ;
            int t = 0;
            int sign = 1;
            for (int x=0;x<image.Width;x+=res)
            {
                int y0 = 0;
                if(sign == -1)
                {
                    y0 = image.Height-1;
                }
                for (int y = y0; (y < image.Height)&&(sign==1) || (y >= 0) && (sign == -1); y+=sign*res)
                {
                    Color pix = image.GetPixel(x, y);
                    r = 2.0 * (pix.R/255.0) - 1.0;
                    v = 2.0 * (pix.G / 255.0) - 1.0;
                    b = 2.0 * (pix.B / 255.0) - 1.0;
                    px = 2.0*(x/(double)image.Width)-1.0;
                    py = 2.0*(y / (double)image.Height) - 1.0;
                    for(int k=0;k<nbRepet;k++)
                    {
                        double r_=r;
                        double v_ = v;
                        double b_ = b;
                        double px_ = px;
                        double py_ = py;
                        double pr = (double)k / (double)nbRepet;
                        if(t>=1)
                        {
                            r_ = retour[0, t - 1] * (1.0 - pr) + r * pr;
                            v_ = retour[1, t - 1] * (1.0 - pr) + v * pr;
                            b_ = retour[2, t - 1] * (1.0 - pr) + b * pr;
                            px_ = retour[3, t - 1] * (1.0 - pr) + px * pr;
                            py_ = retour[4, t - 1] * (1.0 - pr) + py * pr;
                        }
                        retour[0, t] = r_;
                        retour[1, t] = v_;
                        retour[2, t] = b_;
                        retour[3, t] = px_;
                        retour[4, t] = py_;
                        t++;
                    }
                }
                sign *= -1;
            }
            return retour;
        }
        public static Bitmap ImageEcriture(double[,] signal, int x, int y,double scale)
        {
            Bitmap bp = new Bitmap(x,y);
            int xt = 0;
            int yt = 0;
            for(int i=0;i<signal.GetLength(1);i++)
            {
                int y_ = (int)(Math.Max(0, Math.Min(1, 0.5 * ((signal[4, i] / scale) + 1.0))) * (double)(y-1));
                int x_ = (int)(Math.Max(0, Math.Min(1, 0.5 * ((signal[3, i] / scale) + 1.0))) * (double)(x-1));
                int b = (int)(Math.Max(0, Math.Min(1, 0.5 * ((signal[2, i] / scale) + 1.0))) * (double)255);
                int v = (int)(Math.Max(0, Math.Min(1, 0.5 * ((signal[1, i] / scale) + 1.0))) * (double)255);
                int r = (int)(Math.Max(0, Math.Min(1, 0.5 * ((signal[0, i] / scale) + 1.0))) * (double)255);
                double dist = Math.Sqrt((xt - x_) * (xt - x_) + (yt - y_) * (yt - y_));
                int jmax = 2*Math.Max((int)dist, 1);
                for (int j=1;j<=jmax;j++)
                {
                    double k =(j / (double)jmax);
                    int xr = (int)(k * x_ + (1.0 - k) * xt);
                    int yr = (int)(k * y_ + (1.0 - k) * yt);
                    Color pix = bp.GetPixel(xr, yr);
                    if (Color.Equals(pix,Color.FromArgb(255,0,0,0)))
                    {
                        bp.SetPixel(xr, yr, Color.FromArgb(r, v, b));
                    }
                    else
                    {
                        
                        int r0 = pix.R;
                        int v0 =pix.G;
                        int b0 =pix.B;
                        bp.SetPixel(xr, yr, Color.FromArgb((r+r0)/2, (v+v0)/2, (b+b0)/2));
                    }
                }
                
                xt = x_;
                yt = y_;
            }
            return bp;
        }
        //Signal
        public static int nexpPow2(int N)
        {
            return (int)Math.Pow(2, Math.Floor(Math.Log(N, 2)) + 1);
        }
        public static double sinc(double x)
        {
            //Sinus cardinal de x
            if (x == 0)
            {
                return 1;
            }
            else
            {
                return Math.Sin(x) / x;
            }
        }
        public static double[] DSP(double[] x, int N)
        {
            //Calcul la dsp de x, en N points, entre -fe et fe
            int n = x.GetLength(0);
            double[] retour = new double[N];
            double[] retour_ = new double[N];
            for (int i = 0; i < N; i++)
            {
                if (i % 100 == 0)
                    Console.WriteLine(i + "/" + N);
                double re = 0;
                double im = 0;
                for (int k = 0; k < n; k++)
                {
                    re += x[k] * Math.Cos(-2.0 * Math.PI * (double)k * (double)i / (double)N);
                    im += x[k] * Math.Sin(-2.0 * Math.PI * (double)k * (double)i / (double)N);
                }
                retour[i] = (re * re + im * im) / ((double)N);
            }
            for (int i = N / 2; i < N; i++)
            {
                retour_[i - N / 2] = retour[i];
                retour_[i] = retour[i - N / 2];
            }
            return retour_;
        }
        public static double[] Conv(double[] u, double[] v)
        {
            //Convolue u et v
            int n = u.Length;
            int p = v.Length;
            double[] res = new double[n + p - 1];
            for (int i = 0; i < n + p + -1; i++)
            {
                double s = 0;
                for (int j = 0; j < p; j++)
                {
                    int z = i - p + 1 + j;
                    if (z < n && z >= 0)
                    {
                        s += v[j] * u[z];
                    }
                }
                res[i] = s;
            }
            return res;
        }
        public static double[] CentrerConv(double[] c, int n)
        {
            //Centre le retour d'un produit de convolution pour qu'il ait une taille n
            double[] retour = new double[n];
            int delta = (c.Length - n) / 2;
            for (int i = 0; i < n; i++)
            {
                retour[i] = c[i + delta];
            }
            return retour;
        }
        public static double[] Filtrer(double[] x, int Ordre, double f0, double deltaf, double fe)
        {
            //Filtre le signal x par bande centrée sur f0 et de largeur deltaf, avec un ordre donné, sachant la frequence d'echantillonage fe 
            double[] filtre = new double[Ordre];
            for (int i = 0; i < Ordre; i++)
            {
                double t = ((double)i - (Ordre / 2.0)) / fe;
                double cos = Math.Cos(2.0 * Math.PI * f0 * t);
                double sincard = sinc(Math.PI * deltaf * t) * (deltaf / fe);
                filtre[i] = cos * sincard;
            }
            return CentrerConv(Conv(x, filtre), x.Length);
        }
        public static double[] TransposerFreq(double[] x, double fp, double fe)
        {
            //Transpose le signal x sur fp et -fp, sachant fe la frequence d'echantillonage
            double[] ret = copiedb(x);
            for (int i = 0; i < ret.Length; i++)
            {
                double t = ((double)i - (ret.Length / 2.0)) / fe;
                double cos = Math.Cos(2.0 * Math.PI * fp * t); ;
                ret[i] = ret[i] * cos;
            }
            return ret
;
        }
        public static double[] Transmettre(double[,] messages, double fe)
        {
            //Transforme un ensemble de signaux (lignes) en un seul signal temporel modulé en amplitude, sachant la frequence d'echantillonage fe
            int NbMessages = messages.GetLength(0);
            int LongMessages = messages.GetLength(1);
            double deltaf = fe / (2.0 * NbMessages);
            double[] temp = Vect1Val(0, LongMessages);
            for (int i = 0; i < NbMessages; i++)
            {
                double fp = (i + 0.5) * deltaf;
                temp = sommeVect(temp, TransposerFreq(Filtrer(ligne(messages, i),100,0,deltaf,fe), fp, fe));
            }
            return temp;
        }
        public static double[,] Recevoir(double[] signal, int NbMessages, int ordreFiltre, double fe)
        {
            //Inverse de la fonction Transmettre, signal non lissé
            return Recevoir(signal, NbMessages, ordreFiltre, fe, true, 0);
        }
        public static double[,] Recevoir(double[] signal, int NbMessages, int ordreFiltre, double fe, bool Analogique, double Ts)
        {
            //Inverse de la fonction Transmettre, si non Analogique, le signal est lissé entre {-1, 1} 
            int LongMessages = signal.GetLength(0);
            int p = LongMessages;
            double[,] retour = new double[NbMessages, p];
            double deltaf = fe / (2.0 * NbMessages);
            for (int i = 0; i < NbMessages; i++)
            {
                Console.WriteLine(i);
                double fp = (i + 0.5) * deltaf;
                double[] ligne = Filtrer(signal, ordreFiltre, fp, deltaf, fe);
                ligne = TransposerFreq(ligne, fp, fe);
                ligne = Filtrer(ligne, ordreFiltre, 0, deltaf, fe);
                if (!Analogique)
                {
                    ligne = BinVersAnalog(AnalogVersBin(ligne, Ts, fe), Ts, fe);
                }
                for (int q = 0; q < p; q++)
                {
                    retour[i, q] = ligne[q];
                }
            }
            return retour;
        }
        public static double[] BinVersAnalog(bool[] bits, double Ts, double fe)
        {
            //Convertit un tableau de bits en signal analogique de periode TS, pour une frequence d'echantillonage fe
            int Ns = (int)(Ts * fe);
            int n = Ns * bits.Length;
            double[] retour = new double[n];
            int k = 0;
            for (int i = 0; i < bits.Length; i++)
            {
                double val;
                if (bits[i])
                {
                    val = 1;
                }
                else
                {
                    val = -1;
                }
                for (int j = 0; j < Ns; j++)
                {
                    retour[k] = val;
                    k++;
                }
            }
            return retour;
        }
        public static bool[] AnalogVersBin(double[] sig, double Ts, double fe)
        {
            //Convertit un signal analogique de periode TS en tableau de bool, pour une frequence d'echantillonage fe
            int Ns = (int)(Ts * fe);
            int n = sig.Length / Ns;
            bool[] retour = new bool[n];
            int k = 0;
            for (int i = 0; i < n; i++)
            {
                double val = 0;
                for (int j = 0; j < Ns; j++)
                {
                    val += sig[k];
                    k++;
                }
                retour[i] = val > 0;
            }
            return retour;
        }
        public static double PuissanceSignal(double[] sig, double fe)
        {
            //Calcule la puissance su signal sig, fe = frequence echantillonage
            //Parseval
            return Math.Pow(normeVect(sig), 2) / sig.Length;
        }
        public static double[] BruiterSignal(double[] sig, double RSB, double fe, ref Random r)
        {
            //Bruite le signal avec un bruit additif gaussien pour un certain rsb
            double ps = PuissanceSignal(sig, fe);
            double pb = ps / Math.Pow(10, RSB / 10.0);
            double sigma = Math.Sqrt(pb);
            double[] result = new double[sig.Length];
            for (int i = 0; i < sig.Length; i++)
            {
                double r1 = Math.Max(r.NextDouble(), double.Epsilon);
                double r2 = r.NextDouble();
                double rg = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Sin(2.0 * Math.PI * r2) * sigma;
                result[i] = sig[i] + rg;
            }
            return result;
        }
        private static double[] copiedb(double[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }
        private static double[] sommeVect(double[] V1, double[] V2)
        {
            //renvoie une somme de V1 et V2
            int Taille = Math.Min(V1.Length, V2.Length);
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i] + V2[i];
            }
            return Resultat;
        }
        private static double normeVect(double[] V)
        {
            //Norme du vecteur v
            return Math.Sqrt(produitScalaire(V, V));
        }
        private static double produitScalaire(double[] V1, double[] V2)
        {
            //Produit scalaire de deux vecteurs
            int Taille = Math.Min(V1.Length, V2.Length);
            double Somme = 0.0;
            for (int i = 0; i < Taille; i++)
            {
                Somme += V1[i] * V2[i];
            }
            return Somme;
        }
        public static double[] ligne(double[,] M, int i)
        {
            //Recupere la ligne i de M
            int k = M.GetLength(1);
            double[] ret = new double[k];
            for (int l = 0; l < k; l++)
            {
                ret[l] = M[i, l];
            }
            return ret;
        }
        public static double[] Vect1Val(double val, int n)
        {
            //Applique f a chaque terme de x
            double[] res = new double[n];
            for (int i = 0; i < n; i++)
            {
                res[i] = val;
            }
            return res;
        }
        public static double[,] copieMat(double[,] M)
        {
            int n = M.GetLength(0);
            int k = M.GetLength(1);
            double[,] res = new double[n, k];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    res[i, j] = M[i, j];
                }
            }
            return res;
        }
        public static double[,] vectVersMat(double[] v)
        {
            //transforme un vecteur en matrice
            double[,] retour = new double[v.Length, 1];
            for (int i = 0; i < v.Length; i++)
            {
                retour[i, 0] = v[i];
            }
            return retour;
        }
        public static void enregistrerMatrice(double[,] M, string Path)
        {
            //Enregistre une matrice
            using (System.IO.StreamWriter file =
           new System.IO.StreamWriter(@Path))
            {
                for (int i = 0; i < M.GetLength(0); i++)
                {
                    string s = "";
                    for (int j = 0; j < M.GetLength(1); j++)
                    {
                        s = s + M[i, j].ToString() + ";";
                    }
                    file.WriteLine(s);
                }
            }
        }
    }
}
