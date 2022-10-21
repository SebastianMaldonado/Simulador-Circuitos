/*
|====================================================================|
*                            |Interfaz|         
* Descripción:                                                        
*   Clase utilizada para la gestión de todas las clases encargadas
*   del tratamiento del circuito, de igual manera, se encarga
*   de organizar las entradas y la interfaz para el usuario
|====================================================================|
*/

public class Interfaz {

  private Area_trabajo trabajo;
  private Circuito circuito;
  
  private int menu;
  private int pnt_x, pnt_y;
  
  private ArrayList<Boton> bt_portada;
  private ArrayList<Boton> bt_creacion;
  private ArrayList<Boton> bt_visualizacion;
  private ArrayList<Boton> bt_error;
  
  //Opciones de Visualización
    private boolean corriente_ramal;
    private boolean corriente_malla;
    private boolean circuito_desconexo;
    private boolean nombre_nodos;
    private boolean info_ramal;
    
  
  
  public Interfaz (int tam_x, int tam_y) {
    this.pnt_x = tam_x;
    this.pnt_y = tam_y;
    
    this.trabajo = new Area_trabajo();
    this.circuito = new Circuito();
    this.circuitoPred1 ();
    
    construirPortada();
    construirCreacion();
    construirVisualizacion();
    
    
    this.menu = 2;
    
    this.corriente_ramal = true;
    this.corriente_malla = true;
    this.circuito_desconexo = true;
    this.nombre_nodos = true;
    this.info_ramal = true;
  }
  
  public void inicializarMotor () {
    trabajo.motor.inicializar();
  }
  
  //---------------------------------------------------------------------------\\
  //                        |Circuitos Predeterminados|                        \\
  //---------------------------------------------------------------------------\\  
  
  //---------------------------|Circuito 1 - Cuadrado|---------------------------\\
  private void circuitoPred1 () {
    Nodo nodo1 = new Nodo(1000, 1000);
    Nodo nodo2 = new Nodo(2000, 1000);
    Nodo nodo3 = new Nodo(1000, 2000);
    Nodo nodo4 = new Nodo(2000, 2000);
    
    this.circuito.agregarNodo(nodo1);
    this.circuito.agregarNodo(nodo2);
    this.circuito.agregarNodo(nodo3);
    this.circuito.agregarNodo(nodo4);
    this.circuito.agregarRamal(new Ramal(nodo1, nodo2));
    this.circuito.agregarRamal(new Ramal(nodo2, nodo4));
    this.circuito.agregarRamal(new Ramal(nodo4, nodo3));
    this.circuito.agregarRamal(new Ramal(nodo3, nodo1));
  }
   
  
  //---------------------------------------------------------------------------\\
  //                        |Visualización de Interfaz|                        \\
  //---------------------------------------------------------------------------\\  
  
  //---------------------------|Visualizar Interfaz|---------------------------\\
  public void visualizar (float ratonX, float ratonY) {
    switch (menu) {
      case 1:    //Portada
        visualizarBotonesInterfaz (ratonX, ratonY);
        break;
      case 2:    //Área de Trabajo (Creación)
        trabajo.mostrarAreaTrabajo(circuito);
        trabajo.detectarSeleccion(ratonX, ratonY, circuito);
        trabajo.mostrarVistaPrevia(ratonX, ratonY, circuito);
        trabajo.mostrarModificacionComp ();
        visualizarBotonesInterfaz (ratonX, ratonY);
        break;
      case 3:    //Área de Trabajo (Visualización)
        trabajo.mostrarAreaTrabajo(circuito);
        
        if (this.nombre_nodos)
          trabajo.mostrarNombreNodos(circuito);
        
        if (this.info_ramal)
          trabajo.mostrarInfoRamal(ratonX, ratonY, circuito);
        
        visualizarBotonesInterfaz (ratonX, ratonY);
        break;
    
    }
  }
  
  //---------------------------|Visualizar Botones de Interfaz|---------------------------\\
  private void visualizarBotonesInterfaz (float ratonX, float ratonY) {
     switch (menu) {
      case 1:    //Portada
        visualizarBotones(this.bt_portada, ratonX, ratonY, true);
        break;
      case 2:    //Área de Trabajo (Creación)
        visualizarBotones(this.bt_creacion, ratonX, ratonY, true);
        break;
      case 3:    //Área de Trabajo (Visualización)
        visualizarBotones(this.bt_creacion, ratonX, ratonY, false);
        visualizarBotones(this.bt_visualizacion, ratonX, ratonY, new boolean[] {corriente_ramal, corriente_malla, circuito_desconexo, nombre_nodos, info_ramal});
        break;
    }
  }
  
  //---------------------------|Visualizar Grupo de Botones|---------------------------\\
  private void visualizarBotones (ArrayList<Boton> botones, float ratonX, float ratonY, boolean bloq) {
    for (int i = 0; i < botones.size(); i++) {
      if (!bloq) {  //Si el botón está bloquedado
        botones.get(i).mostrar_bloq();
      } else {
        if (botones.get(i).sobreBoton(ratonX, ratonY))
          botones.get(i).mostrar_selc();
        else
          botones.get(i).mostrar();
      }
    }
  }
  
  private void visualizarBotones (ArrayList<Boton> botones, float ratonX, float ratonY, boolean[] bloq) {
    for (int i = 0; i < botones.size(); i++) {
      if (bloq.length > i)
        if (!bloq[i]) {  //Si el botón está bloquedado
          botones.get(i).mostrar_bloq();
          continue;
        }
        
      if (botones.get(i).sobreBoton(ratonX, ratonY))
        botones.get(i).mostrar_selc();
      else
        botones.get(i).mostrar();
    }
  }
  
  
  //---------------------------------------------------------------------------\\
  //                         |Actividad según Eventos|                         \\
  //---------------------------------------------------------------------------\\  
  
  //---------------------------|Click Izquierdo|---------------------------\\
  public void clickIzquierdo (float ratonX, float ratonY) {
    switch (menu) {
      case 1:    //Portada
        accionarBotonesPortada(ratonX, ratonY);
        break;
      case 2:    //Área de Trabajo (Creación)
        accionarBotonesCreacion(ratonX, ratonY);
        
        if (buscarBotonAccionado(this.bt_creacion, ratonX, ratonY) < 0)  //Si no fue presionado ningún botón
          trabajo.armarCircuito(ratonX, ratonY, circuito);
        break;
      case 3:    //Área de Trabajo (Visualización)
        accionarBotonesVisualizacion(ratonX, ratonY);
        break;    
    }
  }
  
  //---------------------------|Click Derecho|---------------------------\\
  public void clickDerecho (float ratonX, float ratonY) {
    switch (menu) {
      case 1:    //Portada
        break;
      case 2:    //Área de Trabajo (Creación)
        trabajo.cambiarOrientacion();
        break;
      case 3:    //Área de Trabajo (Visualización)
        break;
    
    
    }
  }
  
  //---------------------------|Presionar Tecla|---------------------------\\
  public void presionarTecla (char tecla) {
    switch (menu) {
      case 1:    //Portada
        break;
      case 2:    //Área de Trabajo (Creación)
        trabajo.motor.mover();
        trabajo.seleccionarModo(tecla);
        trabajo.ingresarDatoComponente(tecla);
        break;
      case 3:    //Área de Trabajo (Visualización)
        trabajo.motor.mover();
        break;
    
    }
  }
  
  
  //---------------------------------------------------------------------------\\
  //                     |Actividad según Grupo de Botones|                    \\
  //---------------------------------------------------------------------------\\  
  
  private int buscarBotonAccionado(ArrayList<Boton> botones, float ratonX, float ratonY) {
    if (botones == null)
      return -1;
    
    for (int i = 0; i < botones.size(); i++)
      if (botones.get(i).sobreBoton(ratonX, ratonY))
        return i;
        
    return -1;
  }
  
  private void accionarBotonesPortada (float ratonX, float ratonY) {
    int boton = buscarBotonAccionado(this.bt_portada, ratonX, ratonY);
    
    switch (boton) {
    }
  }
  
  private void accionarBotonesCreacion (float ratonX, float ratonY) {
    int boton = buscarBotonAccionado(this.bt_creacion, ratonX, ratonY);
    
    switch (boton) {
      case 0:
        trabajo.seleccionarModo ('q');
        break;
      case 1:
        trabajo.seleccionarModo ('s');
        break;
      case 2:
        trabajo.seleccionarModo ('c');
        break;
      case 3:
        trabajo.seleccionarModo ('r');
        break;
      case 4:
        trabajo.seleccionarModo ('f');
        break;
      case 5:
        trabajo.seleccionarModo ('e');
        break;
      case 6:
        trabajo.seleccionarModo ('q');
        Resolvedor_circuitos r = new Resolvedor_circuitos();
        r.resolverCircuito(this.circuito);
        menu = 3;
        break;
    }
  }
  
  private void accionarBotonesVisualizacion (float ratonX, float ratonY) {
    int boton = buscarBotonAccionado(this.bt_visualizacion, ratonX, ratonY);
    
    switch (boton) {
      case 0:
        this.corriente_ramal = !this.corriente_ramal;
        break;
      case 1:
        this.corriente_malla = !this.corriente_malla;
        break;
      case 2:
        this.circuito_desconexo = !this.circuito_desconexo;
        break;
      case 3:
        this.nombre_nodos = !this.nombre_nodos;
        break;
      case 4:
        this.info_ramal = !this.info_ramal;
        break;
      case 5:
        this.circuito.limpiarCircuito ();
        menu = 2;
        break;
    }
  }
  
  
  //---------------------------------------------------------------------------\\
  //                         |Ajustes de Dimensiones|                          \\
  //---------------------------------------------------------------------------\\  
  
  //---------------------------|Ajuste de Dimensiones|---------------------------\\
  //dim [true] dimensión x | [false] dimensión y
  public int ajustarDim (int pos, boolean dim) {
    if (dim)
      return (pos * this.pnt_x)/1920;
      
    return (pos * this.pnt_y)/1080;
  }
  
  
  //---------------------------------------------------------------------------\\
  //                        |Construcción de Interfaz|                         \\
  //---------------------------------------------------------------------------\\
  
  //---------------------------|Crear Botones Portada|---------------------------\\
  private void construirPortada() {
    this.bt_portada = new ArrayList<Boton>();
  }
  
  //---------------------------|Crear Botones del Menú de Creación|---------------------------\\
  private void construirCreacion() {
    this.bt_creacion = new ArrayList<Boton>();
    
    int x = ajustarDim(1820, true);
    int tamx = ajustarDim(100, true);
    int tamy = ajustarDim(100, false);
    
    for (int i = 0; i <= 5; i++) {    //Crear Botones de Modos de Creación
      crearBoton(this.bt_creacion, "Modo" + "-" + (i + 1), x, ajustarDim(200 + (100 * i), false), tamx, tamy);
    }
    
    //Botón de Resolución del Circuito
    crearBoton(this.bt_creacion, "Resolver", ajustarDim(1400, true), ajustarDim(850, false), ajustarDim(450, true), ajustarDim(175, false));
  }
  
  //---------------------------|Crear Botones Visualización|---------------------------\\
  private void construirVisualizacion() {
    this.bt_visualizacion = new ArrayList<Boton>();
    
    int y = ajustarDim(925, false);
    int tamx = ajustarDim(100, true);
    int tamy = ajustarDim(100, false);
    
    for (int i = 0; i <= 4; i++) {    //Crear Botones de Visualización
      crearBoton(this.bt_visualizacion, "Ver" + "-" + (i + 1), ajustarDim(50 + (125 * i), true), y, tamx, tamy);
    }
    
    //Botón de Volver a Creación de Circuito
    crearBoton(this.bt_visualizacion, "Volver", ajustarDim(1400, true), ajustarDim(850, false), ajustarDim(450, true), ajustarDim(175, false));
  }
  
  //---------------------------|Crear Botón en Grupo|---------------------------\\
  private void crearBoton (ArrayList<Boton> grupo, String nombre, int x, int y, int tx, int ty) {
    grupo.add(new Boton(x, y, tx, ty));
    grupo.get(grupo.size() - 1).ingresarImagen(loadImage(nombre + ".png"), loadImage(nombre + "-1.png"), loadImage(nombre + "-2.png"));
  }
}
