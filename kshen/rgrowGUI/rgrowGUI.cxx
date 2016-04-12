
//FLTK
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Group.H>
#include <FL/fl_draw.H>
#include <stdio.h>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Scroll.H>
#include <iostream>
// COMPILE: fltk-config --use-images --compile mywindow.cxx


//show cordinates


// Callback: when user picks 'Quit'


class MyDesk : public Fl_Group {
private:
    int xcoord, ycoord, width, height;
    

public:

    MyDesk(int X, int Y, int W, int H, const char *L=0) : Fl_Group(X,Y,W,H,L) {
        xcoord = X;
        ycoord = Y;
        width = W;
        height = H;

        color(48);

    }

    int handle(int e) {
        int ret = Fl_Group::handle(e);
        
        switch ( e ) {
            case FL_ENTER:
                ret = 1;                // FL_ENTER: must return(1) to receive FL_MOVE
                break;
            case FL_MOVE:               // FL_MOVE: mouse movement causes 'user damage' and redraw..
                damage(FL_DAMAGE_USER1);
                ret = 1; 
                break;
            case FL_PUSH:
              char c[30];
              sprintf(c, "Seed location: x=%d y=%d", Fl::event_x(), Fl::event_y());
              fl_message(c);
              
                
        }
       return ret;
        
    }


    // Draw mouse coords in small black rectangle
    void draw_coords() {
        // Coordinates as a string
        char s[80];
        sprintf(s, "x=%d y=%d", Fl::event_x(), Fl::event_y());
        // Black rect
        fl_color(FL_BLACK);
        fl_rectf(10 ,10,200,25);
        // White text
        fl_color(FL_WHITE);
        fl_font(FL_HELVETICA, 18);
        fl_draw(s, 15 , 25);
    }
    void draw() {
        // User damage ONLY? just draw coords and done
        if ( damage() == FL_DAMAGE_USER1 ) {
            draw_coords();
            return;
        }
        
        // Draw coords last
        draw_coords();
    }
    
};






int main(int argc, char *argv[]) {

    if ( argc != 2 ) {
        std::cout<<"usage: "<< argv[0] <<" <filename>\n";
        return 1;
    }
        
    fl_register_images();
    Fl_Window     win2(1024,768, "SegmentationGUI");                     
    Fl_Shared_Image *img = Fl_Shared_Image::get(argv[1]);

   
                    // make a window
    if(!img)
    {
        std::cout<<"Open File Error!\n";
        return 1;
    }


    Fl_Box  box(0,0,img->w(),img->h());     // widget that will contain image
    box.image(img);                             // attach jpg image to box

    MyDesk desk(0,0,img->w(),img->h());
    win2.resizable(win2);
    win2.show();
    
    
    return(Fl::run());
} 
    