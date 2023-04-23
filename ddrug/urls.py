from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path
from rest_framework_simplejwt.views import TokenObtainPairView, TokenRefreshView


from .views import  (DrugListView, DrugCardView,detailDrug, DrugCreateView, DrugUpdateView,#createDrug, updateDrug, 
    detailVitekcard,  VitekcardListView, VitekastListView, Importhandler_VITEK, smartsQuery, 
    ketcher_test,iframe_url, API_VITEK_ASTList, API_Drug_List, MIC_COADDListView, MIC_COADDCardView, 
    MIC_PubListView, MIC_PubCardView, MIC_PubListView, MIC_PubCardView)#VitekcardListView,



urlpatterns = [
    # API path
    path('api/token/', TokenObtainPairView.as_view(), name='token_obtain_pair'),
    path('api/token/refresh/', TokenRefreshView.as_view(), name='token_refresh'),
    path('api-vitek-ast/', API_VITEK_ASTList.as_view(), name="api_vitekast"),
    path('api-drug/', API_Drug_List.as_view(), name="api_drug"),
    # path('vitek-ast/create/', VITEK_ASTCreate.as_view()),
    # path('vitek-ast/<pk>/', VITEK_ASTUpdate.as_view()),
    # path('vitek-ast/<pk>/delete/', VITEK_ASTDelete.as_view()),
    # Normal path
    path('drug_card', DrugCardView.as_view(), name="drug_card"),
    path('drug_list', DrugListView.as_view(), name="drug_list"),
    path('drug/<str:pk>', detailDrug, name="drug_detail"),
    path('drug_detail_structure/<str:pk>', smartsQuery, name="smartsquery"),
    path('createDrug/', DrugCreateView.as_view(), name="drug_create"),
    path('updateDrug/<str:pk>', DrugUpdateView.as_view(), name="drug_update"),
    path('vitekcard_list', VitekcardListView.as_view(), name="vitekcard_list"),
    path('vitekast_list', VitekastListView.as_view(), name="vitekast_list"),
    path('vitekcard_detail/<str:pk>', detailVitekcard, name="vitekcard_detail"),

    path('mic-coadd_list', MIC_COADDListView.as_view(), name="mic_coadd_list"),
    path('mic-coadd_card', MIC_COADDCardView.as_view(), name="mic_coadd_card"),
    path('mic-pub_list', MIC_PubListView.as_view(), name="mic_pub_list"),
    path('mic-pub_card', MIC_PubCardView.as_view(), name="mic_pub_card"),
     

    path("import/<str:process_name>/", Importhandler_VITEK.as_view(), name="import-VITEK"),
    path("ketcher_test/", ketcher_test, name="ketcher_test"),
    path("ketcher/", iframe_url, name="ketcher"),

    
]