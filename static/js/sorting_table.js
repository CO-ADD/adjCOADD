$(document).ready(function () {

    $(".order_field").click(function () {

        var order_name = $(this).text().toString()
        // console.log(order_name);
        if (localStorage.getItem('dorder') !== null) {
            order_name = order_name
            localStorage.removeItem('dorder')
        } else {
            localStorage.setItem('dorder', '-')
            order_name = localStorage.getItem('dorder') + order_name
            console.log(order_name)
        }
        $("#order_input").val(order_name);
        $("#django-filter-siderbar").submit();

    })

});