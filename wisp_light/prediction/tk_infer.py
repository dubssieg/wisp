import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from wisp_light.prediction.predict import TaxoPredictor


class PredictApp:
    def __init__(self, master):
        self.master = master
        master.title("Interface MicroTaxo")

        self.model_dir = tk.StringVar()
        self.params_file = tk.StringVar()
        self.fna_path = tk.StringVar()
        self.predict_dir = tk.StringVar()
        self.raw_pred = tk.BooleanVar(value=False)
        self.predictor = None  # ← Initialisation du predictor (persistent)

        self.build_layout()

    def build_layout(self):
        # Frame gauche : entrées
        left_frame = tk.Frame(self.master)
        left_frame.grid(row=0, column=0, padx=10, pady=10, sticky="n")

        # Frame droite : affichage des résultats
        right_frame = tk.Frame(self.master)
        right_frame.grid(row=0, column=1, padx=10, pady=10, sticky="n")

        row = 0
        tk.Label(left_frame, text="Répertoire du modèle").grid(row=row, column=0, sticky='e')
        tk.Entry(left_frame, textvariable=self.model_dir, width=50).grid(row=row, column=1)
        tk.Button(left_frame, text="Parcourir", command=self.browse_model_dir).grid(row=row, column=2)

        row += 1
        tk.Label(left_frame, text="Fichier de paramètres YAML").grid(row=row, column=0, sticky='e')
        tk.Entry(left_frame, textvariable=self.params_file, width=50).grid(row=row, column=1)
        tk.Button(left_frame, text="Parcourir", command=self.browse_params_file).grid(row=row, column=2)

        row += 1
        tk.Label(left_frame, text="Fichier .fna").grid(row=row, column=0, sticky='e')
        tk.Entry(left_frame, textvariable=self.fna_path, width=50).grid(row=row, column=1)
        tk.Button(left_frame, text="Parcourir", command=self.browse_fna_file).grid(row=row, column=2)

        row += 1
        tk.Label(left_frame, text="Répertoire de prédiction (optionnel)").grid(row=row, column=0, sticky='e')
        tk.Entry(left_frame, textvariable=self.predict_dir, width=50).grid(row=row, column=1)
        tk.Button(left_frame, text="Parcourir", command=self.browse_predict_dir).grid(row=row, column=2)

        row += 1
        tk.Checkbutton(left_frame, text="Activer raw_pred", variable=self.raw_pred).grid(row=row, column=1, sticky='w', pady=5)


        row += 1
        tk.Button(left_frame, text="Charger le modèle", command=self.load_predictor,
                  bg="blue", fg="white").grid(row=row, column=1, pady=5)

        row += 1
        tk.Button(left_frame, text="Lancer la prédiction", command=self.run_prediction,
                  bg="green", fg="white").grid(row=row, column=1, pady=10)

        # Zone de texte pour les résultats
        tk.Label(right_frame, text="Résultats (table)").pack(anchor='w')
        self.result_text = scrolledtext.ScrolledText(
            right_frame, width=140, height=30, wrap=tk.NONE,
            font=("Courier", 10), background="#f9f9f9"
        )
        self.result_text.pack()

    def load_predictor(self):
        model_dir = self.model_dir.get()
        params_file = self.params_file.get()
        predict_dir = self.predict_dir.get() or None

        if not model_dir or not params_file:
            messagebox.showerror("Erreur", "Spécifiez le modèle et le fichier de paramètres.")
            return

        try:
            self.predictor = TaxoPredictor(model_dir, params_file, predict_dir=predict_dir)
            messagebox.showinfo("Succès", "Modèle chargé avec succès.")
        except Exception as e:
            messagebox.showerror("Erreur", f"Erreur lors du chargement du modèle :\n{str(e)}")

    def browse_model_dir(self):
        path = filedialog.askdirectory()
        if path:
            self.model_dir.set(path)

    def browse_params_file(self):
        path = filedialog.askopenfilename(filetypes=[("YAML files", "*.yaml *.yml")])
        if path:
            self.params_file.set(path)

    def browse_fna_file(self):
        path = filedialog.askopenfilename(filetypes=[("FNA files", "*.fna")])
        if path:
            self.fna_path.set(path)

    def browse_predict_dir(self):
        path = filedialog.askdirectory()
        if path:
            self.predict_dir.set(path)

    def run_prediction(self):
        if self.predictor is None:
            messagebox.showwarning("Modèle non chargé", "Veuillez d'abord charger le modèle.")
            return

        fna_path = self.fna_path.get()
        raw_pred = self.raw_pred.get()

        if not fna_path:
            messagebox.showerror("Erreur", "Veuillez spécifier un fichier .fna.")
            return

        try:
            self.result_text.delete('1.0', tk.END)
            self.predictor.predict_one_fna(fna_path, raw_pred=raw_pred, verbose=False)
            table_str = self.predictor.results_table.to_string(index=False)
            self.result_text.insert(tk.END, table_str)
            messagebox.showinfo("Succès", "Prédiction terminée.")
        except Exception as e:
            self.result_text.insert(tk.END, f"Erreur : {str(e)}\n")
            messagebox.showerror("Erreur", str(e))

def main():
    root = tk.Tk()
    app = PredictApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
# if __name__ == "__main__":
#     root = tk.Tk()
#     app = PredictApp(root)
#     root.mainloop()
