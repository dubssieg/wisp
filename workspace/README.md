# Cahier des charges

## Développement (~2 mois)

### Prise en main (~0.5 mois)

- [ ] Script pour créer une installation propre (dossiers requis, plus de simplicité pour la localisation fichier de paramètres)
- [ ] Intégrer le script permettant de faire les graphiques à l'outil en ligne de commande (dossier `scripts`)
- [ ] Intégrer le script permettant de télécharger les génomes de référence à l'outil en ligne de commande (dossier `scripts`)

### Améliorations (~1.5 mois)

- [ ] Mettre les derniers paramètres de modèle dans le fichier de paramètres
- [ ] Trouver un moyen de rejeter moins de lectures (pour le moment, ce qui est en dessous du seuil est rejeté : penser à de l'assemblage en amont ?)
- [ ] Effectuer les prédictions par batches et plus une à une
- [ ] Améliorer le système de rapports pour que ceux-ci soient lisibles lorsque l'on passe à l'échelle (récupération des valeurs dans un fichier log de résumé par exemple)

## Tests (~0.5 mois)

- [ ] Tests unitaires sur les fonctions d'encodage, de gestion des arbres taxonomiques
- [ ] Vérifier les arborescences de dossiers et automatiser leur construction

## Validation (~3.5 mois)

Les résultats issus de la validation peuvent être source de modifications sur des fonctionnalités (si on découvre que la fiabilité de l'outil n'est pas impactée par l'usage de lectures de 6kb au lieu de 10kb, nul besoin d'effectuer un pré-assemblage). Les génomes de la *refseq bacteria* ont d'ores et déjà été téléchargés, me demander l'accès.

- [ ] Validation des paramètres et définition des paramètres par défaut
- [ ] Explorer si utiliser [l'alphabet complet](https://international.neb.com/tools-and-resources/usage-guidelines/single-letter-codes) (implémenté) donne des meilleurs résultats que l'alphabet simple (A, T, C et G) en modifiant les arguments de complémentarité, et en gérant les erreurs possibles associées.
- [ ] Valider que le programme passe à l'échelle (pour le moment, base de données 14k espèces et jeu de données 1 million de lectures passe)
- [ ] Si impossible de rejeter moins de lectures, trouver des paramètres efficaces pour des lectures moyennes (~6kb)
- [ ] Valider les résultats avec des données test (fragments de séquences issues des séquences de référence) puis des jeux de données mock connus, et enfin des vraies données MinION
