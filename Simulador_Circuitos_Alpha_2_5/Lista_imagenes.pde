/*
|====================================================================|
*                        |Lista de Imágenes|
* Descripción:                                                        
*   Lista enlazada circular encargada de almacenar las imágenes 
*   para utilizar en el motor gráfico
*   Ingrese una imagen de tipo .png a partir de su nombre
|====================================================================|
*/

class Lista_imagenes {
  PImage imagen;      //Objeto imagen
  String archivo;     //Nombre archivo
  boolean ajustado;   //[verdadero] Si ha sido ajustado | [falso] si debe serlo
  boolean cargar;     //[verdadero] Si ha sido cargado | [falso] si debe serlo
  
  Lista_imagenes siguiente;
  
  
  //----------------------------|Buscar una imagen|----------------------------//
  Lista_imagenes encontrar (String nombre) {
    Lista_imagenes temp = this;
    boolean encontrado = false;
    
    while ((temp.siguiente != this) && (!encontrado)){
      temp = temp.siguiente;
      
      if (temp.archivo == nombre)
        encontrado = true;
    }
    
    if (encontrado) {  //Si se encontró
      return temp;       //Enviar imagen
    } else {           //Si NO se encontró
      ingresar(nombre);    //Crear imagen
      return temp.siguiente;       
    }
  }
  
  
  //----------------------------|Ingresar una imagen|----------------------------//
  //Crea un nuevo objeto en la lista e ingresa los datos correspondientes
  void ingresar (String archivo) {
    if (this.imagen == null) {  //Ingreso de imagen por defecto
      this.archivo = archivo;
      this.ajustado = false;
      this.siguiente = this;
      
      this.imagen = this.cargar();
    } else {                    //Ingreso del resto de imágenes
      Lista_imagenes temp = this;
      
      while (temp.siguiente != this) {
        temp = temp.siguiente;
      }
      
      Lista_imagenes nuevo = new Lista_imagenes();
      temp.siguiente = nuevo;
      
      nuevo.archivo = archivo;
      nuevo.ajustado = false;
      nuevo.siguiente = this;
      
      nuevo.imagen = nuevo.cargar();
      
      temp = null;
      System.gc();
    }
  }
  
  
  //----------------------------|Cargar una imagen|----------------------------//
  //Carga la imagen desde el directorio del juego
  PImage cargar () {
    PImage temp;
    temp = loadImage (this.archivo + ".png");
    
    if (temp == null)
      temp = loadImage ("defecto.png");
    
    return temp;
  }
  
  
  //----------------------------|Ajustar una imagen|----------------------------//
  //Calcula la proporción de los objetos en pantalla
  void ajustar (float prop, float tm_x, float tm_y) {
    if (!ajustado) {  //Si no ha sido ajustada
      //Generar Tamaño Relativo
        int tx = (int)(tm_x * prop);
        int ty = (int)(tm_y * prop);
        
      //Si la imagen se agranda más del 25% de la proporción actual
        if ((tm_x < tx * 1.5) && (tm_y < ty * 1.5))
          this.imagen = this.cargar();  //Volver a cargar imagen
        
      this.imagen.resize(tx, ty);
      this.ajustado = true;
    }
  }
  
  
  //----------------------------|Permitir el ajuste|----------------------------//
  //Cambia la variable ajustado de todas las imágenes de la lista
  //Por lo que cuando se muestren serán ajustadas
  void reajustar () {
    Lista_imagenes temp = this;
    temp.ajustado = false;
    
    while (temp.siguiente != this){
      temp = temp.siguiente;
      temp.ajustado = false;
    }
  }
  
  
  //----------------------------|Volver a cargar las Imágenes|----------------------------//
  //Vuelve a cargar las imágenes para ajustar la resolución
  void refrescar () {
    Lista_imagenes temp = this;
    
    while (temp.siguiente != this){
      temp = temp.siguiente;
      int x = temp.imagen.width;
      int y = temp.imagen.height;
      
      temp.imagen = temp.cargar();
      
      temp.imagen.resize(x, y);
    }
    
    temp = null;
    System.gc();
  }
}
