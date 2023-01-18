
  $(document).ready(function () {
    $(document).on("dblclick", ".editablerow", function () {
        console.log("editable...")
       var value=$(this).text();
       var input="<input type='text' class='input-data' value='"+value+"' />";
       $(this).html(input);
       $(this).removeClass("editablerow");
    });
   
    $(document).on("blur", ".input-data", function(){
         var value=$(this).val();
         var td= $(this).parent("td");
         $(this).remove();
         td.html(value);
         td.addClass("editablerow");
    });
    $(document).on("keypress", ".input-data", function(e){
     
      if(e.keyCode===13 ){

        var value=$(this).val();
        var td= $(this).parent("td");
        $(this).remove();
        td.html(value);
        td.addClass("editablerow");
      }
    })

    $(document).on("click", ".sendToserver", function () {
        console.log("to server...")

        var $row = $(this).closest("tr"),       // Finds the closest row <tr> 
            $td_2 = $row.find("td:nth-child(2)").data("value").toString();
            $td_3= $row.find("td:nth-child(3)").text().toString();
            $td_4= $row.find("td:nth-child(4)").text().toString();
            data={"dict_value": $td_2, "dict_class":$td_3, "dict_desc":$td_4}
            console.log(data)
            sendToServer(data)
    });

    const csrftoken = getCookie('csrftoken');
    const sendToServer = (data) => {
    $.ajax({
        url: '/dict_update', 
        type: "POST",
        headers: { 'X-CSRFToken': csrftoken },
        data: data,

    })
        .done((response) => {
            console.log(response)
        
        })
        .fail(() => {
            console.log("Error occured")
        })
}

  })

