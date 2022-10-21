  /*
|-------------------------------//-------------------------------|
|  Simulador de Circuitos - Alpha 2.4   (20/10/2022)             |
|  Proyecto Independiente                                        |
|                                                                |
|  Desarrolladores:                                              |
|                  - Yordi Gonzales                              |
|                  - Juan David Maestre                          |
|                  - Sebastián Maldonado                         |
|                  - Jaime Vergara                               |
|-------------------------------//-------------------------------|
 
  Descripción de versión:
      
*/

//---------------------------------------------------------------------------\\
//                            |Variables Globales|                           \\
//---------------------------------------------------------------------------\\
import java.util.*;
import org.apache.commons.math3.linear.*;
Interfaz interfaz;


void setup() {
  fullScreen();
  interfaz =  new Interfaz(width, height);
  interfaz.inicializarMotor();
}


void draw() {
  interfaz.visualizar(mouseX, mouseY);
}

//---------------------------------------------------------------------------\\
//                            |Eventos del Ratón|                            \\
//---------------------------------------------------------------------------\\
void mousePressed() {
  if (mouseButton == LEFT)
    interfaz.clickIzquierdo(mouseX, mouseY);
  
  if (mouseButton == RIGHT)
    interfaz.clickDerecho(mouseX, mouseY);
}

//---------------------------------------------------------------------------\\
//                           |Eventos del Teclado|                           \\
//---------------------------------------------------------------------------\\
void keyPressed() {
  interfaz.presionarTecla(key);
}
