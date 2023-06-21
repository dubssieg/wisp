# Cahier des charges

## Développement

- [ ] Script pour créer une installation propre (dossiers requis, plus de simplicité pour la localisation fichier de paramètres)
- [ ] Intégrer le script permettant de faire les graphiques à l'outil en ligne de commande (dossier `scripts`)
- [ ] Intégrer le script permettant de télécharger les génomes de référence à l'outil en ligne de commande (dossier `scripts`)
- [ ] Mettre les derniers paramètres de modèle dans le fichier de paramètres
- [ ] Trouver un moyen de rejeter moins de lectures (pour le moment, ce qui est en dessous du seuil est rejeté : penser à de l'assemblage en amont ?)
- [ ] Effectuer les prédictions par batches et plus une à une

## Tests

- [ ] Tests unitaires sur les fonctions d'encodage, de gestion des arbres taxonomiques
- [ ] Vérifier les arborescences de dossiers et automatiser leur construction

## Validation

- [ ] Validation des paramètres et définition des paramètres par défaut
- [ ] Explorer si utiliser [l'alphabet complet](https://international.neb.com/tools-and-resources/usage-guidelines/single-letter-codes) (implémenté) donne des meilleurs résultats que l'alphabet simple (A, T, C et G)
- [ ] Valider que le programme passe à l'échelle (pour le moment, base de données 14k espèces et jeu de données 1 million de lectures passe)
- [ ] Si impossible de rejeter moins de lectures, trouver des paramètres efficaces pour des lectures moyennes (~6kb)
- [ ] Valider avec des données test (fragments de séquences issues des séquences de référence) puis des jeux de données mock connus, et enfin des vraies données MinION
