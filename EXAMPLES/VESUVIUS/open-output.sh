if [[ "$OSTYPE" == "linux-gnu" ]]; then
        echo 'linux-gnu'
        xdg-open concplot_cum.pdf &
        xdg-open concplot.pdf &
        xdg-open gsd.pdf &
        xdg-open parxplot.pdf &
elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo 'Mac OSX'
        open concplot_cum.pdf  
        open concplot.pdf  
        open gsd.pdf  
        open parxplot.pdf
elif [[ "$OSTYPE" == "cygwin" ]]; then
        echo 'cygwin'
elif [[ "$OSTYPE" == "msys" ]]; then
        echo 'msys'
elif [[ "$OSTYPE" == "win32" ]]; then
        echo 'win32'
elif [[ "$OSTYPE" == "freebsd"* ]]; then
        echo 'freebsd'
else
        echo 'unknown'
fi
