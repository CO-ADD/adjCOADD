function load_Indicator(trigger, initial, loader, initial_switcher = 'on'){
    // if (initial_switcher === 'off'){
    //     initial.addClass("not-visible")
    // }
    if (trigger.is('form')){
        trigger.submit(function(){     
            initial.addClass("not-visible")
            loader.removeClass("not-visible")})

    }else{

        trigger.click(function(){     
            initial.addClass("not-visible")
            loader.removeClass("not-visible")
            
        }) 
    }

}