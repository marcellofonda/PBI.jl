# Radiografia in contrasto di fase

I programmi in questo repository rappresentano un po' di strumenti che ho sviluppato per la mia tesi sull'imaging ai raggi X in contrasto di fase. Il focus della tesi è sul contrasto a propagazione libera.

## Le basi fisiche della trattazione

Le assunzioni e relative approssimazioni principali che vengono usate nella trattazione sono:

1. "Projection approximation" che assume che il campione sia "sottile", cioè che tutti i suoi effetti sulla radiazione possano essere accumulati su un piano ortogonale alla propagazione della radiazione (il piano del campione, o "object plane"). Il campione può quindi essere descritto come una funzione a due variabili, ossia tramite una densità planare di "attività ottica". [^1] Questa funzione avrà due componenti: una di attenuazione/assorbimento e una di traslazione fase, che chiamiamo $\Phi$.
2. "Paraxial approximation", che assume che i vettori d'onda che descrivono la radiazione vengano deviati, durante l'interazione con il campione, di un angolo molto piccolo, cosicché l'immagine registrata possa essere considerata come una funzione locale della densità di attività ottica: ogni pixel sull'immagine dipende solo da un piccolo intorno del corrispondente punto sul piano del campione.
3. "Weak phase contrast": gli effetti dovuti all'interzaione sono deboli. [^dubbio] Questo permette di approssimare al primo ordine in $\Phi$ gli effetti di contrasto di fase.
4. "Monomorficity condition": il campione è composto da un unico materiale, di cui sono noti $\delta$ e $\beta$ ($n=\delta+i\beta$ è l'indice di rifrazione complesso del materiale). Ciò permette di scrivere $\Phi(x,y)=-k\delta t(x,y)$ ed il termine di attenuazione come $e^{-2k\beta t(x,y)}$, con $t(x,y)$ lo spessore dell'oggetto nel piano del campione e $k$ il numero d'onda della radiazione.

L'equazione che risulta per l'intensità da queste approssimazioni è

$$
I(x,y)=I_0(1-\frac{z_1\delta}{2k\beta}\nabla^2)e^{-2k\beta t(x,y)},
$$

dove $z_1$ è la distanza tra oggetto e rilevatore.


[^1]: Mi sono inventato io il termine "attività ottica" in questo contesto, non so quanto sia adeguato.

[^dubbio]: Non ho capito perché in realtà, siccome sfruttiamo proprio il fatto che il laplaciano diverga per avere gli effetti di bordo...


## I programmi

Ci sono tre programmi principali:
* `integratore.jl` è l'espressione della "projection approximation": dato un modello tridimensionale in formato `.obj`, ne calcola lo spessore in ogni punto del piano dell'oggetto ( $t(x,y)$ ) e salva il risultato come un'immagine in bianco e nero. La proiezione viene fatta, per scelta personale, sul piano $yz$ invece che sul piano $xy$, perché di solito i raggi si propagano orizzontalmente ed è più comodo modellare oggetti sul piano orizzontale. Questa scelta, per quanto significativa nella produzione del modello, è ininfluente per tutto il resto.

	Le tre funzioni principali sono:
	* `load_obj(filename::AbstractString)`

	Data una `AbstractString` contenente il nome del file, ne salva i vertici in un array (`vertices`) di vettori a tre coordinate e le facce in un array (`faces`) di vettori a tre indici, relativi ciascuno al corrispettivo vertice del triangolo in `vertices`. Si assume infatti che ogni faccia sia un triangolo. Restituisce i due array `vertices` e `faces`.
	* `project_mesh(vertices, faces, n)`

	Dati gli array `vertices` e `faces`, salva, per ogni punto di una griglia `n`$\times$`n`, le coordinate $x$ di tutti i punti della superficie della mesh che hanno coordinate $y$ e $z$ corrispondenti al punto della griglia. La griglia ha dimensioni scelte in modo da contenere perfettamente il modello. Se non ci sono vertici corrispondenti al punto della griglia, viene salvato il punto della superficie del triangolo la cui proiezione contiene il punto della griglia. Se il punto appartiene a un triangolo ortogonale al piano del campione (la sua proiezione è degenere), viene salvata due volte la coordinata $x$ del baricentro del triangolo. In questo modo si può tener conto dell'effettiva esistenza di tale intersezione, ma si possono anche evitare problemi nel calcolo dello spessore del modello (un numero dispari di punti registrati per ogni pixel potrebbe essere interpretato come una mesh non chiusa. Dovendo integrare una funzione sulla mesh, inoltre, si avrebbero problemi con i domini).

	Restituisce una griglia di vettori, ciascuno contenente un numero pari di coordinate $x$, in ordine crescente
	* `thickness(vector)`

	Dato un vettore di numeri, ne calcola la somma a segni alterni. Se questi numeri rappresentano gli estremi di segmenti, il risultato è la lunghezza totale dei segmenti. Applicato al risultato di `project_mesh`, restituisce la mappa di spessori $t(x,y)$ del modello 3D importato.

* `radiografo.jl` è l'espressione della formula finale per l'intensità. Data un'immagine rappresentante $t(x,y)$, calcola la radiografia risultante, su un sensore avente la stessa risoluzione dell'immagine in input. Ciò crea un po' di confusione, siccome per un'immagine decente al contrasto di fase è richiesta una certa dimensione per i pixel del rilevatore (dell'ordine del $\mu$m), le dimensioni del campione risultano direttamente proporzionali alla risoluzione dell'immagine che lo rappresenta. Una soluzione a questo problema sarà più immediata nel momento in cui implementerò una convoluzione che renda indipendenti la dimensione del pixel del rilevatore, la risoluzione dell'immagine e le dimensioni del campione.

* `filtrolaplaciano.jl` è un file dimostrativo dell'applicazione del laplaciano discreto. Legge un'immagine e la salva in bianco e nero, affiancata al suo laplaciano. Niente di più e niente di meno.




## Note importanti

I file `.obj` devono contenere delle mesh che siano:
* _triangolate_: tutte le facce devono essere triangoli. Non sono ammessi poligoni
* _chiuse_: la superficie della mesh non deve avere bordo.


## TODO
* `integratore.jl`: Aggiungere un controllo sul fatto che la mesh sia triangolata: dovrebbe dare errore se incontra una faccia non triangolare.
* Aggiungere un sistema per gestire le dimensioni fisiche dei pixel in modo consapevole: per ora c'è una correlazione troppo automatica e variabile con le dimensioni del modello e la risoluzione della proiezione. Bisognerebbe renderli indipendenti e gestire
* Stabilizzare la tipizzazione. Questo potrebbe rivelarsi utile in caso di implementazione su scheda grafica, ma anche migliorare le performance in generale.
* Creare un programma unico che legga da `.obj` e restituisca le radiografie pronte. In altre parole, utilizzare i programmi attuali come pacchetti, lasciando al loro interno soltanto le definizioni delle funzioni.
* Implementare un algoritmo di convoluzione che permetta di gestire risoluzioni diverse senza impazzire. Nello specifico, di poter caratterizzare il rilevatore con una funzione che, tramite convoluzione,  possa dare radiografie accurate anche con risoluzioni sempre maggiori nell'output di `integratore.jl`.
