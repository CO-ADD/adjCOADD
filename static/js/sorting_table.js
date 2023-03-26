$(document).ready(function () {
    if (sessionStorage.getItem('dorder') !== null) {
        $(".order_field").toggleClass("desc")
    }
    else {
        $(".order_field").toggleClass("asc")
    }
    $(".order_field").dblclick(function () {

        var order_name = $(this).text().toString()
        if (sessionStorage.getItem('dorder') !== null) {
            order_name = order_name
            sessionStorage.removeItem('dorder')
            // $(this).toggleClass(".asc")
            $(this).toggleClass(".desc")
        } else {
            sessionStorage.setItem('dorder', '-')
            order_name = sessionStorage.getItem('dorder') + order_name
            console.log(order_name)
        }
        $("#order_input").val(order_name);
        $("#django-filter-siderbar").submit();

    })

});