import tkinter as tk
import sympy as sp

class Matrix:
    def __init__(self, name=None, matrix=sp.Matrix.eye(3)):
        self.name = name
        self.matrix = matrix
        self.text = self.convert()
    def setMatrix(self, matrix):
        self.matrix = matrix
        self.text = self.convert()
    def convert(self):
        stringVal = str(self.matrix)
        open = False
        mat = ''
        for s in stringVal:
            if s == '[':
                open = True
                continue
            elif s == ']':
                open = False
                mat += '\n'
                continue
            elif s == ',':
                continue
            elif open:
                mat += s
        return mat.strip()
    def __str__(self):
        return self.name

class LinAlg:
    def __init__(self, master):
        master.title("Linear Algebra")
        matrices = tk.Frame(master)
        matrices.pack(side=tk.LEFT)

        mats = [Matrix("A"), Matrix("B"), Matrix("C")]
        choice = tk.StringVar(matrices)
        choice.set(mats[0])

        self.curMatrix = mats[0]
        menu = tk.OptionMenu(matrices, choice, *mats, command=self.chooseMatrix)
        menu.grid(row=0)
        self.txt = tk.Text(matrices, width=30, height=15)
        self.txt.insert("1.0", mats[0].convert())
        self.txt.bind("<Leave>", self.parseMatrix)
        self.txt.grid(row=1)

        btns = tk.Frame(master)
        btns.pack(side=tk.LEFT)

        btnTranspose = tk.Button(btns, text="X^t", command=self.transpose)
        btnTranspose.grid(row=0, sticky="ew")
        btnInverse = tk.Button(btns, text="X^-1", command=self.inverse)
        btnInverse.grid(row=1, sticky="ew")
        btnREF = tk.Button(btns, text="REF(X)", command=self.ref)
        btnREF.grid(row=2, sticky="ew")
        btnRank = tk.Button(btns, text="Rank(X)", command=self.rank)
        btnRank.grid(row=3, sticky="ew")
        btnDet = tk.Button(btns, text="Det(X)", command=self.det)
        btnDet.grid(row=4, sticky="ew")
        btnLU = tk.Button(btns, text="LU(X)", command=self.lu)
        btnLU.grid(row=5, sticky="ew")
        btnCol = tk.Button(btns, text="C(X)", command=self.col)
        btnCol.grid(row=0, column=1, sticky="ew")
        btnRow = tk.Button(btns, text="R(X)", command=self.row)
        btnRow.grid(row=1, column=1, sticky="ew")
        btnNull = tk.Button(btns, text="N(X)", command=self.null)
        btnNull.grid(row=2, column=1, sticky="ew")
        btnDiag = tk.Button(btns, text="Diag(X)", command=self.diag)
        btnDiag.grid(row=3, column=1, sticky="ew")
        btnEig = tk.Button(btns, text="Eig(X)", command=self.eig)
        btnEig.grid(row=4, column=1, sticky="ew")
        btnCofact = tk.Button(btns, text="Cofact(X)", command=self.cofact)
        btnCofact.grid(row=5, column=1, sticky="ew")

        sol = tk.Frame(master)
        sol.pack(side=tk.LEFT)

        solLbl = tk.Label(sol, text="Solution:")
        solLbl.grid(row=0)
        self.solTxt = tk.Text(sol, width=30, height=15)
        self.solTxt.grid(row=1)

    def transpose(self):
        trans = Matrix(matrix=self.curMatrix.matrix.transpose())
        self.showSolution(trans.text)

    def inverse(self):
        try:
            inv = Matrix(matrix=self.curMatrix.matrix.inv())
            self.showSolution(inv.text)
        except:
            self.showSolution("not invertible")
    
    def ref(self):
        ref = Matrix(matrix=self.curMatrix.matrix.echelon_form())
        self.showSolution(ref.text)

    def det(self):
        try:
            det = self.curMatrix.matrix.det()
            self.showSolution(str(det))
        except:
            self.showSolution("not square")

    def rank(self):
        rank = self.curMatrix.matrix.rank()
        self.showSolution(str(rank))

    def col(self):
        self.showSolution(self.showSpace(self.curMatrix.matrix.columnspace()))

    def row(self):
        row = self.curMatrix.matrix.rowspace()
        rt = []
        for r in row:
            rt.append(sp.Matrix.transpose(r))
        self.showSolution(self.showSpace(rt))

    def null(self):
        null = self.curMatrix.matrix.nullspace()
        self.showSolution(self.showSpace(null))
    
    def lu(self):
        L, U, _ = self.curMatrix.matrix.LUdecomposition()
        l = Matrix(matrix=L)
        u = Matrix(matrix=U)
        self.showSolution("L\n" + l.text + "\n\nU\n" + u.text)

    def diag(self):
        try:
            P, D = self.curMatrix.matrix.diagonalize()
            p = Matrix(matrix=P)
            p_inv = Matrix(matrix=P.inv())
            d = Matrix(matrix=D)
            self.showSolution("P\n" + p.text + "\n\nP^-1\n" + p_inv.text + "\n\nD\n" + d.text)
        except:
            self.showSolution("not square")
    
    def eig(self):
        try:
            eig = self.curMatrix.matrix.eigenvects()
            val = ""
            for e in eig:
                val += "Eigenvalue: " + str(e[0])
                if e[1] > 1:
                    val += "\tM=" + str(e[1]) + "\n"
                    val += self.showSpace(e[2]) + "\n\n"
                    continue
                val += '\n'
                vect = Matrix(matrix=e[2])
                val += vect.text[8:] + "\n\n"
            self.showSolution(val)
        except:
            self.showSolution("not square")

    def cofact(self):
        c = self.curMatrix.matrix.cofactor_matrix()
        self.showSolution(Matrix(matrix=c).text)

    def showSpace(self, space):
        m = []
        for s in space:
            m.append(sp.Matrix.transpose(s))
        m = sp.BlockMatrix(m)
        m = m.transpose()
        return Matrix(matrix=m).text

    def showSolution(self, text):
        self.solTxt.delete("1.0", tk.END)
        self.solTxt.insert("1.0", text)

    def chooseMatrix(self, item):
        self.curMatrix = item
        self.txt.delete("1.0", tk.END)
        self.txt.insert("1.0", self.curMatrix.text)

    def parseMatrix(self, event):
        text = self.txt.get("1.0", "end-1c")
        rows = text.split('\n')
        A = []
        for row in rows:
            A.append([int(n.strip()) for n in row.split()])
        self.curMatrix.setMatrix(sp.Matrix(A))

root = tk.Tk()
la = LinAlg(root)
root.mainloop()